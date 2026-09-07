"""
simulator.py – Functions related to virtual Parsivel sampling of model output.

Model data is represented as ``xr.Dataset`` objects whose variable and
coordinate names have been normalised to the canonical vocabulary defined in
:mod:`pyPIPS.model_config`.  The primary entry point for loading model data is
:func:`load_model_dataset`, which calls :func:`model_config.open_and_normalise`.

No "model_dict" / "grid_dict" / "dis_dict" catch-all dictionaries are used in
the public API.  All functions receive explicit, typed arguments.
"""

from __future__ import annotations

import inspect
import os
import glob
from datetime import datetime
from typing import Sequence

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import xarray as xr
from joblib import Parallel, delayed
from metpy.plots import ctables
from scipy.stats import gamma, uniform
from shapely.geometry import LineString, MultiLineString
from numpy.random import poisson

from . import DSDlib as dsd
from . import PIPS as pips
from . import parsivel_params as pp
from . import plotmodule as pm
from . import thermolib as thermo
from .legacy import disdrometer_module as dis
from . import radarmodule as radar
from .model_config import ModelConfig, get_grid_arrays, open_and_normalize
import pyPIPS.dualpara as dualpol

try:
    import pyart
except ImportError:  # pragma: no cover
    pyart = None

# Optional pyCRMtools imports used only by the ARPS helper functions at the
# bottom of this module.
try:
    from pyCRMtools.modules import utils as CRMutils
    from pyCRMtools.pycaps import arps_read
    _HAS_PYCRMTOOLS = True
except ImportError:   # pragma: no cover
    _HAS_PYCRMTOOLS = False
    CRMutils = None
    arps_read = None


# ---------------------------------------------------------------------------
# Module-level constants derived from Parsivel parameters
# ---------------------------------------------------------------------------

rhoacst = 1.0       # kg m^-3 (approximate sea-level air density)
rhorcst = 1000.0    # kg m^-3 (liquid-water density)
#rhohlcst = 850.0    # kg m^-3 (hail density)
cr = rhorcst * np.pi / 6.0
mur = 1.0 / 3.0     # assumed rain shape parameter (gamma-diameter)
muhl = 1.0 / 3.0     # assumed hail shape parameter (gamma-diameter)
mug = 1.0 / 3.0     # assumed graupel shape parameter (gamma-diameter)

sampling_area_default = pp.parsivel_parameters["sensor_area_mm2"] / 1.0e6  # m^2
sampling_width_default = pp.parsivel_parameters["sensor_width_mm"] / 1.0e3  # m
sampling_length_default = pp.parsivel_parameters["sensor_length_mm"] / 1.0e3  # m

D = pp.parsivel_parameters["avg_diameter_bins_mm"] / 1.0e3  # m
Dl = pp.parsivel_parameters["min_diameter_bins_mm"] / 1.0e3  # m
Dr = pp.parsivel_parameters["max_diameter_bins_mm"] / 1.0e3  # m
Dedges = np.append(Dl, Dr[-1])
bin_width = Dr - Dl


# ===========================================================================
# DSD / sampling utilities
# ===========================================================================

def get_Dmax_index(
    Dr: np.ndarray,
    Dmax: float | None,
) -> tuple[float, int]:
    """
    Return ``(Dmax_m, Dmax_index)`` for the given right-edge diameter array.

    Parameters
    ----------
    Dr : np.ndarray
        Right-edge diameter bin values (m or mm).
    Dmax : float or None
        Maximum diameter in m or mm (but must be consistent with Dr).
        If ``None``, the last bin is used.

    Returns
    -------
    tuple of (Dmax, Dmax_index)
    """
    if Dmax is not None:
        Dmax_index = int(np.searchsorted(Dr, Dmax, side="left"))
    else:
        Dmax_index = int(np.size(Dr)) - 1
        Dmax = float(Dr[Dmax_index])
    return Dmax, Dmax_index


def samplegammaDSD(
    N: float,
    lamda: float,
    alpha: float,
    bins: np.ndarray | None = None,
) -> np.ndarray:
    """
    Randomly sample a gamma DSD given *Nt*, *lamda*, and *alpha*.

    Parameters
    ----------
    N : float
        Total number of particles to sample (not a concentration,
        but an actual number of particles).
    lamda : float
        Slope parameter (m^-1).
    alpha : float
        Shape parameter.
    bins : array-like or None
        If provided, histogram the samples into these diameter bins and
        return ND = counts / bin_width.

    Returns
    -------
    np.ndarray
        Raw diameter samples if *bins* is None, otherwise ND per bin.
    """
    n = int(N)
    if n <= 0 or float(lamda) <= 0.0:
        # No particles to sample or undefined DSD; return empty / zero array.
        if bins is None:
            return np.array([])
        return np.zeros(len(bins) - 1)
    scale = 1.0 / float(lamda)
    shape = float(alpha) + 1.0
    # Hmmm... gamma.rsv expects an integer size, but N is a float.
    # For now, just take int(N) as the number of samples
    s = gamma.rvs(shape, scale=scale, size=n)
    if bins is None:
        return s
    ND_sample, _ = np.histogram(s, bins)
    return ND_sample / (bins[1:] - bins[:-1])


def calc_ND(
    pcount_binned: np.ndarray,
    sampling_volumes_D: np.ndarray,
    Dr: np.ndarray,
    Dl: np.ndarray,
    Dmax: float | None,
) -> np.ndarray:
    """
    Calculate the number density ND (#·m^-3·m^-1) for each diameter bin.

    Parameters
    ----------
    pcount_binned : np.ndarray
        Particle counts per bin.
    sampling_volumes_D : np.ndarray
        Sampling volumes as a function of diameter (m^3).
    Dr, Dl : np.ndarray
        Right- and left-edge bin arrays (m).
    Dmax : float or None
        Maximum diameter (m); used to select bins.

    Returns
    -------
    np.ndarray
    """
    Dmax, Dmax_index = get_Dmax_index(Dr, Dmax)
    return pcount_binned / (
        sampling_volumes_D * (Dr[: Dmax_index + 1] - Dl[: Dmax_index + 1])
    )


def calc_sampling_volumes_D(
    Vt: np.ndarray,
    Dr: np.ndarray,
    Dmax: float | None,
    sampling_interval: float,
    sampling_area: float,
) -> np.ndarray:
    """
    Calculate sampling volumes (m^3) as a function of terminal velocity.

    Parameters
    ----------
    Vt : np.ndarray
        Terminal velocity per diameter bin (m s^-1).
    Dr : np.ndarray
        Right-edge diameter bin array (m).
    Dmax : float or None
        Maximum diameter (m).
    sampling_interval : float
        Accumulation interval (s).
    sampling_area : float
        Sensor cross-sectional area (m^2).

    Returns
    -------
    np.ndarray
    """
    Dmax, Dmax_index = get_Dmax_index(Dr, Dmax)
    return Vt[: Dmax_index + 1] * sampling_interval * sampling_area


def create_random_gamma_DSD(
    Nt: float,
    lamda: float,
    alpha: float,
    Vt: np.ndarray,
    sampling_length: float,
    sampling_width: float,
    Dl: np.ndarray,
    Dmid: np.ndarray,
    Dr: np.ndarray,
    Dmin: float = 0.0,
    Dmax: float | None = None,
    sampling_interval: float = 10.0,
    remove_margins: bool = False,
    verbose: bool = False,
    rhocorrect: bool = False,
    rho: float | None = None,
    mask_lowest: bool = True,
    perturb_vel: bool = True,
    hail: bool = False,
    rhohl: float | None = None,
) -> xr.Dataset:
    """
    Given *Nt*, *lamda*, *alpha* create a spatial DSD sample within a volume.

    Returns an ``xr.Dataset`` with variables ``positions``, ``diameters``,
    ``velocities``, ``ND``, ``margin_mask``, ``pcount_binned``,
    ``sampling_volumes_D``.
    """
    Dmax, Dmax_index = get_Dmax_index(Dr, Dmax)
    if verbose:
        print("Dmax_index =", Dmax_index)

    Vtmax = Vt[Dmax_index]
    sampling_area = sampling_length * sampling_width
    sampling_height = Vtmax * sampling_interval
    sampling_volume = sampling_area * sampling_height
    sampling_volumes_D = calc_sampling_volumes_D(
        Vt, Dr, Dmax, sampling_interval, sampling_area
    )

    if verbose:
        print("sampling height =", sampling_height)
        print("sampling volume =", sampling_volume)

    # From Anna James, use Poisson distribution to determine number of particles in the sampling
    # volume, to avoid problems with small (< 1) expected numbers of particles in the sampling
    # volume, which would always be clipped to zero with just int()
    # n = int(Nt * sampling_volume)
    n = int(poisson(Nt * sampling_volume))

    if verbose:
        print("number concentration =", Nt)
        print("number of particles in sampling volume =", n)

    # Check that n > 0 before proceeding with sampling; if not, return None
    if n <= 0:
        return None

    xpos = uniform.rvs(0.0, sampling_length, (n, 1))
    ypos = uniform.rvs(0.0, sampling_width, (n, 1))
    zpos = uniform.rvs(0.0, sampling_height, (n, 1))

    diameters = samplegammaDSD(n, lamda, alpha)

    if verbose and diameters.size:
        print("min/max diameter in sample =", diameters.min(), diameters.max())
        print("maximum allowed diameter =", Dmax)

    diameter_mask = diameters <= Dmax
    if mask_lowest:
        low_mask = diameters > Dr[1]
        if verbose:
            print("number of particles above the lowest two bins =", low_mask.sum())
        diameter_mask = diameter_mask & low_mask

    # Check that there are still particles left after applying the diameter mask; if not, return None
    if diameter_mask.sum() <= 0:
        if verbose:
            print("No particles left after applying diameter mask; returning None.")
        return None

    diameters = diameters[diameter_mask]
    xpos = xpos[diameter_mask]
    ypos = ypos[diameter_mask]
    zpos = zpos[diameter_mask]

    if verbose:
        print("number within allowable diameter range =", diameter_mask.sum())
        print(
            "min/max particle diameter in truncated sample =",
            diameters.min(), diameters.max(),
        )
    if not hail:
        velocities = pips.calc_empirical_fallspeed(
            diameters * 1000.0, correct_rho=rhocorrect, rho=rho
        )
    else:
        if verbose:
            print("Calculating hail particle velocities with rhohl =", rhohl)
        velocities = pips.calc_empirical_fallspeed_hail(
            diameters, rhohl, correct_rho=rhocorrect, rho=rho)
    if verbose:
        print("min/max particle velocity in sample =", velocities.min(), velocities.max())
    depths = velocities * sampling_interval
    keepers = np.where(zpos.squeeze() - depths <= 0.0)

    xpos = xpos[keepers]
    ypos = ypos[keepers]
    zpos = zpos[keepers]
    diameters = diameters[keepers]
    velocities = velocities[keepers]

    margins_xl = xpos.squeeze() - diameters / 2.0 < 0.0
    margins_xr = xpos.squeeze() + diameters / 2.0 > sampling_length
    margins_yl = ypos.squeeze() - diameters / 2.0 < 0.0
    margins_yr = ypos.squeeze() + diameters / 2.0 > sampling_width
    margin_mask = margins_xl | margins_xr | margins_yl | margins_yr

    if verbose:
        print("number of particles that fall through sampling area =", xpos.size)
        print("number that are margin fallers =", margin_mask.sum())

    if remove_margins:
        if verbose:
            print("Removing margin fallers!")
        xpos = xpos[~margin_mask]
        ypos = ypos[~margin_mask]
        zpos = zpos[~margin_mask]
        diameters = diameters[~margin_mask]
        velocities = velocities[~margin_mask]

    _Dedges = np.append(Dl, Dr[-1])
    pcount_binned, _ = np.histogram(diameters, _Dedges)
    ND = calc_ND(pcount_binned, sampling_volumes_D, Dr, Dl, Dmax)

    positions = np.hstack((xpos, ypos, zpos))

    # Check shapes of arrays before returning
    if verbose:
        print("positions shape:", positions.shape)
        print("diameters shape:", diameters.shape)
        print("velocities shape:", velocities.shape)
        print("margin_mask shape:", margin_mask.shape)
        print("ND shape:", ND.shape)
        print("pcount_binned shape:", pcount_binned.shape)
        print("sampling_volumes_D shape:", sampling_volumes_D.shape)

    # Not returning margin_mask for now because of issues with alignment of dimensions, but may
    # revisit later.
    return xr.Dataset(
        {
            "positions": xr.DataArray(positions, dims=["particle", "coord"]),
            "diameters": xr.DataArray(diameters, dims=["particle"]),
            "velocities": xr.DataArray(velocities, dims=["particle"]),
            # "margin_mask": xr.DataArray(margin_mask, dims=["particle"]),
            "ND": xr.DataArray(ND, dims=["diameter_bin"]),
            "pcount_binned": xr.DataArray(pcount_binned, dims=["diameter_bin"]),
            "sampling_volumes_D": xr.DataArray(sampling_volumes_D, dims=["diameter_bin"]),
        }
    )


def uniquetol(a: np.ndarray, tol: float = 1.0e-3) -> np.ndarray:
    """
    Return a copy of *a* with near-duplicates removed.

    Requires *a* to already be sorted in either increasing or decreasing order.
    """
    b = np.array(a)
    d = np.append(True, np.abs(np.diff(b)))
    return b[d > tol]


def combine_sample_and_model_times(
    model_times: Sequence,
    probe_times: Sequence,
) -> np.ndarray:
    """
    Merge model output times and probe sample times into a sorted unique array.

    Parameters
    ----------
    model_times, probe_times : sequences
        Both should contain comparable objects (floats or datetimes).

    Returns
    -------
    np.ndarray
        Sorted, unique combined array.
    """
    combined = sorted(set(model_times) | set(probe_times))
    return np.array(combined)


# ===========================================================================
# Model dataset helpers
# ===========================================================================

def load_model_dataset(cfg: ModelConfig) -> xr.Dataset:
    """
    Load and normalise a model dataset using *cfg*.

    This is the **primary entry point** for ingesting model data into the
    simulator pipeline.  The returned Dataset has canonical variable/coord
    names and ``dx``/``dy`` stored in ``ds.attrs``.

    Parameters
    ----------
    cfg : ModelConfig
        Built with :func:`model_config.make_config` or directly.

    Returns
    -------
    xr.Dataset
    """
    return open_and_normalize(cfg)


def get_fields_at_time(
    ds: xr.Dataset,
    time_sec: float,
    field_names: Sequence[str],
    level: int = 0,
    time_coord: str = "time",
) -> xr.Dataset:
    """
    Extract a 2-D (y, x) Dataset from *ds* at the model time closest to
    *time_sec*.

    Parameters
    ----------
    ds : xr.Dataset
        Normalised model dataset.
    time_sec : float
        Requested time in seconds (must use the same numeric convention as
        the values stored in ``ds[time_coord]``).
    field_names : sequence of str
        Canonical field names to extract.
    level : int
        Vertical level index (0 = lowest model level).
    time_coord : str
        Name of the time coordinate in *ds*.

    Returns
    -------
    xr.Dataset with one DataArray per requested field.  Names absent from
    *ds* are silently skipped.
    """
    tvals = np.asarray(ds[time_coord].values, dtype=float)
    tidx = int(np.argmin(np.abs(tvals - time_sec)))
    rec = ds.isel({time_coord: tidx})

    out: dict[str, xr.DataArray] = {}
    for name in field_names:
        if name not in rec:
            continue
        da = rec[name]
        if da.ndim >= 3:
            da = da.isel({da.dims[0]: level})
        out[name] = da
    return xr.Dataset(out)


def compute_derived_fields(
    fields: xr.Dataset,
    requested: Sequence[str],
) -> xr.Dataset:
    """
    Add derived fields to a field Dataset.

    Computed fields
    ---------------
    UC   : U wind interpolated to cell centres
    VC   : V wind interpolated to cell centres
    PTE  : equivalent potential temperature (requires P, TH, QV)
    rhoa : air density (requires P, TH, QV)

    Parameters
    ----------
    fields : xr.Dataset
        Dataset of 2-D DataArrays (already extracted at a single time).
    requested : sequence of str
        Names of derived fields to add.  Unrecognised names are ignored.

    Returns
    -------
    The *fields* Dataset with derived variables added where possible.
    """
    new_vars: dict = {}
    if "UC" in requested and "U" in fields:
        U = fields["U"].values
        new_vars["UC"] = 0.5 * (U[:-1, :-1] + U[:-1, 1:])
    if "VC" in requested and "V" in fields:
        V = fields["V"].values
        new_vars["VC"] = 0.5 * (V[:-1, :-1] + V[1:, :-1])
    if "PTE" in requested and {"P", "TH", "QV"}.issubset(fields):
        try:
            new_vars["PTE"] = thermo.calpte(
                fields["P"].values, fields["TH"].values, fields["QV"].values
            )
        except Exception:
            print("Warning: cannot calculate PTE – skipping.")
    if "rhoa" in requested and "rhoa" not in fields:
        if {"P", "TH", "QV"}.issubset(fields):
            try:
                new_vars["rhoa"] = thermo.calrho(
                    fields["P"].values, fields["TH"].values, fields["QV"].values
                )
            except Exception:
                print("Warning: cannot calculate rhoa – skipping.")
    if new_vars:
        fields = fields.assign(new_vars)
    return fields


def read_fields_at_model_times(
    ds: xr.Dataset,
    model_times: Sequence[float],
    scalar_vars: Sequence[str],
    vector_vars: Sequence[str] = (),
    derived_vars: Sequence[str] = (),
    level: int = 0,
    time_coord: str = "time",
) -> list[xr.Dataset]:
    """
    Extract and return a list of field Datasets, one per model time.

    Parameters
    ----------
    ds : xr.Dataset
        Normalised model dataset.
    model_times : sequence of float
        Times (s) at which to sample.
    scalar_vars : sequence of str
        Canonical names of scalar fields.
    vector_vars : sequence of str
        Canonical names of (staggered) vector wind fields.
    derived_vars : sequence of str
        Names of derived quantities to compute via
        :func:`compute_derived_fields`.
    level : int
        Vertical level index.
    time_coord : str

    Returns
    -------
    list of xr.Dataset, one per entry in *model_times*.
    """
    all_vars = list(dict.fromkeys(list(scalar_vars) + list(vector_vars)))
    result = []
    for t_req in model_times:
        fields = get_fields_at_time(
            ds, float(t_req), all_vars, level=level, time_coord=time_coord
        )
        fields = compute_derived_fields(fields, derived_vars)
        result.append(fields)
    return result


# ===========================================================================
# Radar sweep reading
# ===========================================================================

def read_sweeps(
    radardir: str,
    radname: str,
    radstarttimestamp: str,
    radstoptimestamp: str,
    fieldnames: Sequence[str],
    el_req: float,
) -> tuple[list, list]:
    """
    Read CFRadial sweeps and return ``(radarsweeplist, sweeptimelist)``.

    Parameters
    ----------
    radardir : str
        Directory containing the radar netCDF files.
    radname : str
        Radar name substring used in the glob pattern.
    radstarttimestamp : str
        Start time as ``'%Y%m%d%H%M%S'``.
    radstoptimestamp : str
        Stop time as ``'%Y%m%d%H%M%S'``.
    fieldnames : sequence of str
        Field names to read.
    el_req : float
        Requested elevation angle (degrees).

    Returns
    -------
    (radarsweeplist, sweeptimelist) sorted by ascending time.
    """
    radpathlist = glob.glob(os.path.join(radardir, f"*{radname}*nc"))
    start_dt = datetime.strptime(radstarttimestamp, "%Y%m%d%H%M%S")
    stop_dt = datetime.strptime(radstoptimestamp, "%Y%m%d%H%M%S")

    radarsweeplist = []
    sweeptimelist = []

    for radpath in radpathlist:
        sweeptime = radar._getsweeptime(radpath)
        if start_dt <= sweeptime <= stop_dt:
            radarsweep = radar.readCFRadial_pyART(
                el_req, radpath, sweeptime, fieldnames, compute_kdp=False
            )
            radarsweeplist.append(radarsweep)
            sweeptimelist.append(sweeptime)

    sorted_sweeptimelist = sorted(sweeptimelist)
    sorted_radarsweeplist = [
        x
        for _, x in sorted(
            zip(sweeptimelist, radarsweeplist), key=lambda pair: pair[0]
        )
    ]
    return sorted_radarsweeplist, sorted_sweeptimelist


def compute_storm_motion(
    feature_start_time: datetime,
    feature_end_time: datetime,
    feature_start_loc: tuple[float, float],
    feature_end_loc: tuple[float, float],
) -> tuple[float, float]:
    """
    Compute storm-motion vector (m s^-1) from start/end positions.

    Parameters
    ----------
    feature_start_time, feature_end_time : datetime
    feature_start_loc, feature_end_loc : (x_km, y_km)
        Positions in km relative to the radar origin.

    Returns
    -------
    (ustorm, vstorm) in m s^-1.
    """
    deltat = (feature_end_time - feature_start_time).total_seconds()
    ustorm = (feature_end_loc[0] - feature_start_loc[0]) * 1000.0 / deltat
    vstorm = (feature_end_loc[1] - feature_start_loc[1]) * 1000.0 / deltat
    return ustorm, vstorm


# ===========================================================================
# Probe location helpers
# ===========================================================================

def get_probe_geolocs_from_files(
    dis_names: Sequence[str],
    dis_types: Sequence[str],
    dis_filepaths: Sequence[str],
    conv_filepaths: Sequence[str | None],
    starttimes: Sequence,
    stoptimes: Sequence,
) -> list[tuple[float, float]]:
    """
    Return a list of ``(lat, lon)`` tuples, one per probe.

    Currently only PIPS data in netCDF format is supported. PIPS location
    is read from the probe's netCDF dataset attributes.

    Parameters
    ----------
    dis_names : sequence of str
        Probe names (used in error messages).
    dis_types : sequence of str
        One of ``'PIPS'`` (others not yet supported).
    dis_filepaths : sequence of str
        Paths to netCDF disdrometer data files.
    conv_filepaths : sequence of str or None
        Paths to netCDF conventional data files (required for PIPS).
    starttimes, stoptimes : sequences
        Time-window limits (currently unused; included for API compatibility).

    Returns
    -------
    list of (lat, lon) tuples.

    Raises
    ------
    NotImplementedError
        If dis_type is not 'PIPS'.
    """
    geo_locs = []
    for dis_name, dis_type, filepath, conv_filepath in zip(
        dis_names, dis_types, dis_filepaths, conv_filepaths
    ):
        if dis_type == "PIPS":
            try:
                parsivel_ds = xr.load_dataset(filepath)
                # Prefer dataset attributes (already set by the processing
                # pipeline).  The canonical attribute is "location", a string
                # of the form "(lat, lon, alt)" that we parse with eval().
                if "location" in parsivel_ds.attrs:
                    loc_tuple = eval(parsivel_ds.attrs["location"])  # noqa: S307
                    lat, lon = float(loc_tuple[0]), float(loc_tuple[1])
                elif "GPS_lat" in parsivel_ds and "GPS_lon" in parsivel_ds:
                    lat = float(np.nanmean(parsivel_ds["GPS_lat"].values))
                    lon = float(np.nanmean(parsivel_ds["GPS_lon"].values))
                else:
                    raise KeyError(
                        f"Cannot find lat/lon in {dis_name} dataset {filepath}"
                    )
                geo_locs.append((lat, lon))
            except (KeyError, IndexError, ValueError) as e:
                raise ValueError(
                    f"Failed to read lat/lon for PIPS probe {dis_name!r} "
                    f"from {filepath}: {e}"
                )
        else:
            raise NotImplementedError(
                f"Probe type {dis_type!r} is not currently supported. "
                f"Only 'PIPS' data in netCDF format is supported. "
                f"(attempted for probe {dis_name!r})"
            )
    return geo_locs


def get_probe_locs_relative_to_radar(
    probe_geolocs: Sequence[tuple[float, float]],
    radar_sweep,
) -> list[tuple[float, float]]:
    """
    Convert geographic probe locations to Cartesian offsets from the radar.

    Parameters
    ----------
    probe_geolocs : sequence of (lat, lon) tuples
    radar_sweep : py-ART radar object
        Used to determine the radar's lat/lon.

    Returns
    -------
    List of ``(dx_m, dy_m)`` tuples in metres.
    """
    if pyart is None:
        raise ImportError("pyart is required for get_probe_locs_relative_to_radar().")

    rlat = radar_sweep.latitude["data"][0]
    rlon = radar_sweep.longitude["data"][0]

    locs = []
    for lat, lon in probe_geolocs:
        dradx, drady = pyart.core.geographic_to_cartesian_aeqd(lon, lat, rlon, rlat)
        dx = float(dradx[0])
        dy = float(drady[0])
        print(f"  probe offset: dx={dx:.1f} m  dy={dy:.1f} m")
        locs.append((dx, dy))
    return locs


def compute_probe_model_offsets(
    probe_radar_offsets: Sequence[tuple[float, float]],
    reference_radar_offset: tuple[float, float],
    reference_model_offset: tuple[float, float],
) -> list[tuple[float, float]]:
    """
    Transform probe offsets from radar-relative to model-relative coordinates.

    Uses a reference feature (e.g., hook echo tip) to align the radar and model
    coordinate systems.  Computes the model-relative position of each probe
    by shifting the radar-relative positions so that the reference feature
    aligns between radar and model frames.

    Parameters
    ----------
    probe_radar_offsets : sequence of (dx, dy) tuples
        Probe positions relative to the radar in metres.
    reference_radar_offset : (dx, dy) tuple
        Position of a reference feature (e.g., hook echo tip) in radar-relative
        coordinates (metres).
    reference_model_offset : (dx, dy) tuple
        Position of the same reference feature in model-relative coordinates
        (metres).

    Returns
    -------
    List of (x_m, y_m) tuples representing model-relative probe positions
    (metres).
    """
    ref_x_radar, ref_y_radar = reference_radar_offset
    ref_x_model, ref_y_model = reference_model_offset

    offsets = []
    for x_radar, y_radar in probe_radar_offsets:
        x_model = x_radar - ref_x_radar + ref_x_model
        y_model = y_radar - ref_y_radar + ref_y_model
        offsets.append((x_model, y_model))
    return offsets


def get_probe_locs_arps_grid(
    ds: xr.Dataset,
    probe_geo_locs: Sequence[tuple[float, float]],
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    """
    Return probe locations within an ARPS model grid with a map projection.

    Only valid for real-case ARPS runs where the Dataset carries a basemap
    object in ``ds.attrs['bgmap']``.

    Parameters
    ----------
    ds : xr.Dataset
        Normalised ARPS dataset.  Must have ``bgmap``, ``dx``, ``dy``,
        ``xs``, ``ys`` attributes / coordinates.
    probe_geo_locs : sequence of (lat, lon) tuples

    Returns
    -------
    (modloc_list, coord_list)
        Model x/y positions and fractional i/j indices for each probe.
    """
    bgmap = ds.attrs["bgmap"]
    dx = ds.attrs["dx"]
    dy = ds.attrs["dy"]
    xc = np.asarray(ds.coords.get("xs", ds.coords["xc"]).values).ravel()
    yc = np.asarray(ds.coords.get("ys", ds.coords["yc"]).values).ravel()

    modloc_list = []
    coord_list = []
    for lat, lon in probe_geo_locs:
        xloc, yloc = bgmap(lon, lat)
        modloc_list.append(np.array([xloc, yloc]))
        iloc = (xloc - xc[1]) / dx
        jloc = (yloc - yc[1]) / dy
        coord_list.append(np.array([iloc, jloc]))
    return modloc_list, coord_list


# ===========================================================================
# Probe time-series reading and selection
# ===========================================================================

def load_pips_probe_datasets(
    dis_types: Sequence[str],
    dis_filepaths: Sequence[str],
    starttimes: Sequence | None = None,
    stoptimes: Sequence | None = None,
) -> list[xr.Dataset]:
    """
    Load PIPS parsivel_combined datasets from netCDF files.

    The parsivel_combined files already contain both DSD data and conventional
    data resampled to the DSD interval with proper variable names
    (e.g. ``windspdavgvec``, ``winddiravgvec``, ``slowtemp``, ``rho``, etc.).
    No separate conventional file or additional resampling is needed.

    Parameters
    ----------
    dis_types : sequence of str
        Instrument types.  Only ``'PIPS'`` is currently supported.
    dis_filepaths : sequence of str
        Paths to parsivel_combined netCDF files, one per probe.
    starttimes, stoptimes : sequences or None
        Optional per-probe time-window bounds.  Pass ``None`` (the default)
        to load the full time range of each file.

    Returns
    -------
    list of xr.Dataset
        One Dataset per probe.

    Raises
    ------
    NotImplementedError
        If any entry in *dis_types* is not ``'PIPS'``.
    """
    n = len(dis_filepaths)
    if starttimes is None:
        starttimes = [None] * n
    if stoptimes is None:
        stoptimes = [None] * n

    datasets = []
    for dis_type, filepath, starttime, stoptime in zip(
        dis_types, dis_filepaths, starttimes, stoptimes
    ):
        if dis_type != "PIPS":
            raise NotImplementedError(
                f"Probe type {dis_type!r} is not currently supported. "
                f"Only 'PIPS' data in netCDF format is supported."
            )
        ds = xr.load_dataset(filepath)
        if starttime is not None or stoptime is not None:
            ds = ds.sel(time=slice(starttime, stoptime))
        datasets.append(ds)
    return datasets


def select_pips_at_target_times(
    probe_datasets: Sequence[xr.Dataset],
    target_times: Sequence[datetime],
    tolerance_sec: float = 10.0,
) -> list[xr.Dataset]:
    """
    Re-index each PIPS probe Dataset to a set of target times.

    For each probe, the nearest available time is selected for every target
    time.  Target times with no data within *tolerance_sec* are filled with
    NaN (for float variables) or NaT (for datetime variables), so callers
    can detect undeployed periods by checking for NaN in any scalar field
    (e.g. ``ds['slowtemp'].isnull()``).

    Parameters
    ----------
    probe_datasets : sequence of xr.Dataset
        Loaded PIPS parsivel_combined Datasets, e.g. from
        :func:`load_pips_probe_datasets`.
    target_times : sequence of datetime
        Times at which to evaluate probe data (e.g. radar sweep times or
        model output times).
    tolerance_sec : float
        Maximum gap (s) between a target time and the nearest PIPS sample.
        Larger gaps result in NaN fill rather than extrapolation.

    Returns
    -------
    list of xr.Dataset
        One Dataset per probe, re-indexed to *target_times*.
    """
    target_idx = pd.DatetimeIndex(target_times)
    tolerance = pd.Timedelta(seconds=tolerance_sec)
    return [
        ds.reindex(time=target_idx, method="nearest", tolerance=tolerance)
        for ds in probe_datasets
    ]


# ===========================================================================
# Composite bookkeeping helpers
# ===========================================================================

def init_composite(
    ds: xr.Dataset,
    composite_width: tuple[float, float],
    searchbox_width: tuple[float, float],
    tracking_level: float,
    tracking_varname: str,
    tracking_extremum: str,
    gridlims: tuple[float, float, float, float],
) -> xr.Dataset:
    """
    Initialise composite bookkeeping from a normalised model Dataset.

    Parameters
    ----------
    ds : xr.Dataset
        Normalised model dataset (output of :func:`load_model_dataset`).
    composite_width : (width_x_m, width_y_m)
        Full width of the composite box in metres.
    searchbox_width : (width_x_m, width_y_m)
        Full width of the feature-tracking search box in metres.
    tracking_level : float
        Height AGL (m) at which to track the feature.
    tracking_varname : str
        Variable to track (e.g. ``'w'``, ``'vortz'``).
    tracking_extremum : str
        ``'max'`` or ``'min'``.
    gridlims : (xmin, xmax, ymin, ymax)
        Spatial limits within which to search (same units as ds coords).

    Returns
    -------
    xr.Dataset
        Composite metadata consumed by :func:`build_composite`.  Grid
        coordinate arrays are stored as Dataset coordinates; scalar
        parameters are stored in ``attrs``.
    """
    grid = get_grid_arrays(ds)
    xc1d = grid["xc1d"]
    yc1d = grid["yc1d"]
    zc1d = grid["zc1d"]
    xe1d = grid["xe1d"]
    ye1d = grid["ye1d"]
    ze1d = grid["ze1d"]
    dx = grid["dx"]
    dy = grid["dy"]
    compositewidthx, compositewidthy = composite_width
    searchboxwidthx, searchboxwidthy = searchbox_width

    igbgn = int(np.searchsorted(xc1d, gridlims[0], side="left"))
    igend = int(np.searchsorted(xc1d, gridlims[1], side="right"))
    jgbgn = int(np.searchsorted(yc1d, gridlims[2], side="left"))
    jgend = int(np.searchsorted(yc1d, gridlims[3], side="right"))
    gridlimindices = [igbgn, igend + 1, jgbgn, jgend + 1]

    xc1dg = xc1d[igbgn:igend]
    yc1dg = yc1d[jgbgn:jgend]
    xe1dg = xe1d[igbgn:igend]
    ye1dg = ye1d[jgbgn:jgend]
    xckm = xc1dg / 1000.0
    yckm = yc1dg / 1000.0

    print(
        f"Starting/ending grid coordinates: "
        f"igbgn={igbgn:d}, igend={igend:d}, "
        f"jgbgn={jgbgn:d}, jgend={jgend:d}"
    )

    zeagl1d = ze1d - ze1d[0]
    zcagl1d = zc1d - zc1d[0]
    print(zeagl1d, zcagl1d)

    ichw = int(compositewidthx / (2.0 * dx))
    jchw = int(compositewidthy / (2.0 * dy))
    ishw = int(searchboxwidthx / (2.0 * dx))
    jshw = int(searchboxwidthy / (2.0 * dy))

    xckm_comp = np.arange(-ichw, ichw + 1) * dx / 1000.0
    yckm_comp = np.arange(-jchw, jchw + 1) * dy / 1000.0

    if tracking_varname == "w":
        tracking_klvl = int(np.where(zeagl1d >= tracking_level)[0][0])
    else:
        tracking_klvl = int(np.where(zcagl1d >= tracking_level)[0][0])

    print(
        f"Tracking {tracking_extremum} {tracking_varname} "
        f"at height {tracking_level:.1f} m (k={tracking_klvl:d})"
    )

    return xr.Dataset(
        coords={
            "xe_g": ("xe_g", xe1dg),
            "ye_g": ("ye_g", ye1dg),
            "xc_g": ("xc_g", xc1dg),
            "yc_g": ("yc_g", yc1dg),
            "xc_km": ("xc_g", xckm),
            "yc_km": ("yc_g", yckm),
            "ze_agl": ("ze_agl", zeagl1d),
            "zc_agl": ("zc_agl", zcagl1d),
            "xc_comp": ("xc_comp", xckm_comp),
            "yc_comp": ("yc_comp", yckm_comp),
        },
        attrs={
            "igbgn": igbgn,
            "igend": igend + 1,
            "jgbgn": jgbgn,
            "jgend": jgend + 1,
            "ichw": ichw,
            "jchw": jchw,
            "ishw": ishw,
            "jshw": jshw,
            "tracking_klvl": tracking_klvl,
            "tracking_varname": tracking_varname,
            "tracking_extremum": tracking_extremum,
            "tracking_level": tracking_level,
            "dx": dx,
            "dy": dy,
        },
    )


def trackfeature(
    var: np.ndarray,
    searchboxlims: tuple[int, int, int, int] | None = None,
    guesscoords: tuple[int, int] | None = None,
    extremum: str = "max",
    debug: bool = False,
) -> tuple[float, int, int]:
    """
    Locate the extremum of a 2-D field within an optional search box.

    Parameters
    ----------
    var : np.ndarray
        2-D field array (shape ``[ny, nx]``).
    searchboxlims : (ibgn, iend, jbgn, jend) or None
        Index limits of the search subdomain.  If ``None`` the whole array
        is searched.
    guesscoords : (i, j) or None
        Previous feature location used to offset *ibgn*, *jbgn* when
        converting relative to absolute indices.
    extremum : {'max', 'min'}
    debug : bool

    Returns
    -------
    (var_extremum, iref, jref) : (float, int, int)
    """
    if searchboxlims is None:
        ibgn, iend, jbgn, jend = 0, var.shape[1] - 1, 0, var.shape[0] - 1
    else:
        ibgn, iend, jbgn, jend = searchboxlims

    if debug:
        fig, ax = plt.subplots()

    varsearch = var[jbgn: jend + 1, ibgn: iend + 1]
    print(varsearch.shape)

    if debug:
        ax.imshow(varsearch, origin="lower")

    if extremum == "max":
        var_extremum = float(varsearch.max())
        flatindex = int(np.argmax(varsearch))
    else:
        var_extremum = float(varsearch.min())
        flatindex = int(np.argmin(varsearch))

    jrel, irel = np.unravel_index(flatindex, varsearch.shape)

    if debug:
        print(f"irel, jrel = {irel}, {jrel}")
        ax.plot([irel], [jrel], marker="o", color="red")

    if guesscoords is None:
        iref, jref = int(irel), int(jrel)
    else:
        iref = ibgn + int(irel)
        jref = jbgn + int(jrel)

    if debug:
        print(f"iref, jref = {iref}, {jref}")

    return var_extremum, iref, jref


def get_composite_grid(ds: xr.Dataset, comp_info: xr.Dataset) -> xr.Dataset:
    """
    Build a storm-centred composite grid as an ``xr.Dataset``.

    Parameters
    ----------
    ds : xr.Dataset
        Normalised model dataset (provides ``dx``, ``dy``, ``zc``).
    comp_info : xr.Dataset
        Output of :func:`init_composite`.

    Returns
    -------
    xr.Dataset
        Dataset with coordinates ``xc``, ``yc``, ``zc``, ``xe``, ``ye``
        centred on the composite origin, plus ``dx``/``dy`` in attrs.
    """
    ichw = comp_info.attrs["ichw"]
    jchw = comp_info.attrs["jchw"]
    dx = comp_info.attrs["dx"]
    dy = comp_info.attrs["dy"]
    grid = get_grid_arrays(ds)
    zc1d = grid["zc1d"]

    xc_comp = np.arange(-ichw, ichw + 1) * dx
    yc_comp = np.arange(-jchw, jchw + 1) * dy
    xe_comp = np.arange(-ichw - 0.5, ichw + 1.5) * dx
    ye_comp = np.arange(-jchw - 0.5, jchw + 1.5) * dy

    return xr.Dataset(
        coords={
            "xc": xc_comp,
            "yc": yc_comp,
            "zc": zc1d,
            "xe": xe_comp,
            "ye": ye_comp,
        },
        attrs={"dx": dx, "dy": dy},
    )


# ===========================================================================
# Transect intersection – plotting helper (separated from logic)
# ===========================================================================

def _plot_transect_location(
    Zmodplot: np.ndarray,
    xcplot: np.ndarray,
    ycplot: np.ndarray,
    xcorplot: np.ndarray,
    ycorplot: np.ndarray,
    xlocs: np.ndarray,
    ylocs: np.ndarray,
    t: int,
    dxmod: float,
    dymod: float,
    plocs: tuple,
    all_times: np.ndarray,
    model_times_rel: np.ndarray,
    tloc: int,
) -> None:
    """Plot a single transect sample-location panel (side-effect only)."""
    fig = None
    ax = None
    ptype = 2
    xlim = [xlocs.min() - 5000.0, xlocs.max() + 10000.0]
    ylim = [ylocs.min() - 10000.0, ylocs.max() + 10000.0]
    clevels = np.arange(0.0, 85.0, 5.0)
    norm, cmap = ctables.registry.get_with_steps("NWSReflectivity", 0.0, 5.0)

    fig, ax = pm.plotsingle(
        fig, ax, ptype, xcplot, ycplot, xcorplot, ycorplot,
        xlim, ylim, Zmodplot,
        clevels, cmap, norm, clevels, "Z (dBZ)",
        None, False, None, 0, None, None, None, None, None, [2000.0, 2000.0],
    )

    print(f"x, y = {xlocs[t]:.1f}, {ylocs[t]:.1f}")
    ax.plot(xlocs[t], ylocs[t], "ko", ms=4)
    ax.plot(dxmod, dymod, "bo", ms=4)
    ax.plot(xcplot[plocs], ycplot[plocs], "kx", ms=5, alpha=0.5)

    if all_times[t] == model_times_rel[tloc]:
        ax.set_title(f"Time = {all_times[t]:06.0f} (model time)")
    else:
        ax.set_title(f"Time = {all_times[t]:06.0f} (intermediate)")


# ===========================================================================
# Transect grid intersection
# ===========================================================================

def find_transect_grid_intersections(
    ds: xr.Dataset,
    probe_times: Sequence,
    probe_xy_offset: tuple[float, float],
    reference_time: datetime,
    ustorm: float,
    vstorm: float,
    debug: bool = False,
) -> xr.Dataset:
    """
    Find intersections of a single probe transect with the model grid.

    Parameters
    ----------
    ds : xr.Dataset
        Normalized model dataset.  Grid coordinates are read from ``ds.coords``
        via :func:`model_config.get_grid_arrays`.
    probe_times : sequence of datetime
        PSD sample times for the probe.
    probe_xy_offset : (x_m, y_m)
        Model-grid-relative x/y position of the probe at the reference time in meters.
    reference_time : datetime
        Reference time for computing elapsed seconds.
    ustorm, vstorm : float
        Storm-motion components (m s^-1).  Used to advect the grid-relative
        probe position.
    debug : bool

    Returns
    -------
    xr.Dataset with:

    ``x``, ``y``               - positions at all time (including grid intersections),
    ``x_sample``, ``y_sample`` - positions at probe sampling times,

    and coordinates ``all_time`` and ``sample_time``
    (both in seconds from *reference_time*).
    """
    grid = get_grid_arrays(ds)
    xe1d = grid["xe1d"]
    ye1d = grid["ye1d"]

    shapely_grid = MultiLineString(
        [((x, ye1d[0]), (x, ye1d[-1])) for x in xe1d]
        + [((xe1d[0], y), (xe1d[-1], y)) for y in ye1d]
    )

    dxmod, dymod = probe_xy_offset
    sampling_times = np.array(
        [(t - reference_time).total_seconds() for t in probe_times]
    )
    sample_xlocs = dxmod - ustorm * sampling_times
    sample_ylocs = dymod - vstorm * sampling_times

    if debug:
        print("sampling_times =", sampling_times)
        print("sample_xlocs =", sample_xlocs)
        print("sample_ylocs =", sample_ylocs)

    path = LineString(np.c_[sample_xlocs, sample_ylocs])

    difference = path.difference(shapely_grid)
    segments = difference.geoms if hasattr(difference, "geoms") else [difference]

    # Collect (t, x, y) triples from ALL vertices of every sub-segment.
    # Interior vertices of each segment are original sample points that lie
    # within one grid cell; shared endpoints between adjacent segments are the
    # grid-edge crossing points.  Taking all vertices gives the full merged set.
    pts: list[tuple[float, float, float]] = []
    for segment in segments:
        sx, sy = segment.xy
        for xi, yi in zip(sx, sy):
            t = (dxmod - xi) / ustorm
            pts.append((t, xi, yi))

    # Sort by time and remove near-duplicate entries (each grid-edge crossing
    # appears as the last point of one segment and the first of the next).
    pts.sort(key=lambda p: p[0])
    all_times_arr = np.array([p[0] for p in pts])
    xlocs_arr = np.array([p[1] for p in pts])
    ylocs_arr = np.array([p[2] for p in pts])

    keep = np.concatenate(([True], np.abs(np.diff(all_times_arr)) > 1.0e-3))
    all_times = all_times_arr[keep]
    xlocs = xlocs_arr[keep]
    ylocs = ylocs_arr[keep]

    # Snap any all_times entry that is within tolerance of a sampling_time to
    # that exact sampling_time value, so that np.isin comparisons work reliably.
    # Uses broadcasting: diffs shape is (len(sampling_times), len(all_times)).
    diffs = np.abs(all_times[np.newaxis, :] - sampling_times[:, np.newaxis])
    closest_sample_idx = np.argmin(diffs, axis=0)
    min_diffs = diffs[closest_sample_idx, np.arange(len(all_times))]
    snap_mask = min_diffs <= 1.0e-3
    all_times[snap_mask] = sampling_times[closest_sample_idx[snap_mask]]

    if debug:
        print("xlocs =", xlocs)
        print("ylocs =", ylocs)
        print("all_times =", all_times)

    return xr.Dataset(
        {
            "x": xr.DataArray(xlocs, dims=["all_time"]),
            "y": xr.DataArray(ylocs, dims=["all_time"]),
            "x_sample": xr.DataArray(sample_xlocs, dims=["sample_time"]),
            "y_sample": xr.DataArray(sample_ylocs, dims=["sample_time"]),
        },
        coords={
            "all_time": all_times,
            "sample_time": sampling_times,
        },
    )


def get_transect_grid_indices(
    transect_ds: xr.Dataset,
    ds: xr.Dataset,
    reference_time_sec: float,
    model_times: Sequence[float] | None = None,
) -> xr.Dataset:
    """
    Map transect locations to model scalar-point grid indices.

    For every point in the ``all_time`` dimension of *transect_ds*, finds the
    (i, j) index of the model grid cell whose center (scalar point) contains
    that location.  Optionally also computes the index of the model output
    time immediately prior to each transect time.

    Parameters
    ----------
    transect_ds : xr.Dataset
        Output of :func:`find_transect_grid_intersections`.
    ds : xr.Dataset
        Normalized model dataset (provides grid coordinates).
    reference_time_sec : float
        Model time (s) corresponding to the zero of the ``all_time``
        coordinate.  Used to convert transect-relative times to absolute
        model times when *model_times* is given.
    model_times : sequence of float, optional
        Sorted array of model output times (s).  When provided, a
        ``time_idx`` variable is included giving the index of the last
        model output time ≤ each transect time.

    Returns
    -------
    xr.Dataset
        Same ``all_time`` coordinate as *transect_ds*, with variables:

        ``i_idx``    – x (column) index into the model scalar-point grid,
        ``j_idx``    – y (row) index into the model scalar-point grid,
        ``time_idx`` – index into *model_times* (only when supplied).
    """
    grid = get_grid_arrays(ds)
    xe1d = grid["xe1d"]
    ye1d = grid["ye1d"]
    nx = len(xe1d) - 1
    ny = len(ye1d) - 1

    xlocs = transect_ds["x"].values
    ylocs = transect_ds["y"].values
    all_times = transect_ds.coords["all_time"].values

    i_idx = np.clip(np.searchsorted(xe1d, xlocs, side="right") - 1, 0, nx - 1)
    j_idx = np.clip(np.searchsorted(ye1d, ylocs, side="right") - 1, 0, ny - 1)

    data_vars: dict = {
        "i_idx": xr.DataArray(i_idx.astype(int), dims=["all_time"]),
        "j_idx": xr.DataArray(j_idx.astype(int), dims=["all_time"]),
    }

    if model_times is not None:
        mt = np.asarray(model_times, dtype=float)
        abs_times = reference_time_sec + all_times
        time_idx = np.clip(
            np.searchsorted(mt, abs_times, side="right") - 1, 0, len(mt) - 1
        )
        data_vars["time_idx"] = xr.DataArray(time_idx.astype(int), dims=["all_time"])

    return xr.Dataset(data_vars, coords={"all_time": all_times})


def plot_transect_grid_intersections(
    ds: xr.Dataset,
    probe_transect: xr.Dataset,
    grid_idx_ds: xr.Dataset | None = None,
    ax=None,
    title: str | None = None,
):
    """
    Plot the model grid, probe transect path, sample points, and intersection
    points for a single probe.

    Optionally accepts the output of :func:`get_transect_grid_indices` to
    color-code markers and highlight traversed grid cells, making it easy to
    see which probe samples and crossing events belong to the same grid cell.

    Parameters
    ----------
    ds : xr.Dataset
        Normalized model dataset (provides grid coordinates).
    probe_transect : xr.Dataset
        Single-probe output of :func:`find_transect_grid_intersections`.
    grid_idx_ds : xr.Dataset, optional
        Output of :func:`get_transect_grid_indices` for the same probe.
        When provided, markers and traversed-cell borders are color-coded
        by grid cell, cycling through the default matplotlib color cycle.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure/axes is created if not provided.
    title : str, optional
        Axes title.

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    grid = get_grid_arrays(ds)
    xe1d = grid["xe1d"]
    ye1d = grid["ye1d"]

    if ax is None:
        _, ax = plt.subplots()

    # Full grey grid background
    for x in xe1d:
        ax.axvline(x, color="0.8", linewidth=0.5, zorder=1)
    for y in ye1d:
        ax.axhline(y, color="0.8", linewidth=0.5, zorder=1)

    x_sample = probe_transect["x_sample"].values
    y_sample = probe_transect["y_sample"].values

    all_time = probe_transect.coords["all_time"].values
    sample_time = probe_transect.coords["sample_time"].values
    is_crossing = ~np.isin(all_time, sample_time)
    x_int = probe_transect["x"].values[is_crossing]
    y_int = probe_transect["y"].values[is_crossing]

    # Transect path always drawn as a thin blue line
    ax.plot(x_sample, y_sample, "b-", linewidth=1.0, zorder=6, label="transect")

    if grid_idx_ds is None:
        ax.plot(x_sample, y_sample, "bo", ms=4, zorder=6, label="sample")
        ax.plot(x_int, y_int, "rx", ms=6, mew=1.5, zorder=6, label="intersection")
        ax.legend(loc="best", fontsize="small")
    else:
        i_all = grid_idx_ds["i_idx"].values
        j_all = grid_idx_ds["j_idx"].values

        # Build an ordered list of unique (i, j) cells in traversal order
        seen: set = set()
        unique_cells: list[tuple[int, int]] = []
        for i, j in zip(i_all, j_all):
            key = (int(i), int(j))
            if key not in seen:
                seen.add(key)
                unique_cells.append(key)

        palette = [c["color"] for c in plt.rcParams["axes.prop_cycle"]]
        cell_color = {
            cell: palette[k % len(palette)]
            for k, cell in enumerate(unique_cells)
        }

        # Draw a colored border around each traversed cell using the actual
        # local cell-edge coordinates so variable spacing is handled correctly.
        # The full rectangle is drawn at a lower zorder, then the west and
        # south edges are redrawn at a higher zorder so they are always visible.
        for (ci, cj), color in cell_color.items():
            x0, x1 = xe1d[ci], xe1d[ci + 1]
            y0, y1 = ye1d[cj], ye1d[cj + 1]
            rect = Rectangle(
                (x0, y0), x1 - x0, y1 - y0,
                fill=False, edgecolor=color, linewidth=1.5, zorder=2,
            )
            ax.add_patch(rect)
            # West edge (x = x0, y0 → y1) and south edge (y = y0, x0 → x1)
            ax.plot([x0, x0], [y0, y1], color=color, linewidth=1.5, zorder=5)
            ax.plot([x0, x1], [y0, y0], color=color, linewidth=1.5, zorder=5)

        # Grid indices at sample times (sample_time is a subset of all_time)
        at_sample_idx = np.searchsorted(all_time, sample_time)
        i_idx_sample = i_all[at_sample_idx]
        j_idx_sample = j_all[at_sample_idx]

        # Sample points, grouped by cell
        for (ci, cj), color in cell_color.items():
            mask = (i_idx_sample == ci) & (j_idx_sample == cj)
            if mask.any():
                ax.plot(x_sample[mask], y_sample[mask], "o", ms=4,
                        color=color, zorder=6)

        # Grid indices at crossing-only points
        crossing_indices = np.where(is_crossing)[0]
        i_idx_cross = i_all[crossing_indices]
        j_idx_cross = j_all[crossing_indices]

        # Crossing points, grouped by cell
        for (ci, cj), color in cell_color.items():
            mask = (i_idx_cross == ci) & (j_idx_cross == cj)
            if mask.any():
                ax.plot(x_int[mask], y_int[mask], "x", ms=6, mew=1.5,
                        color=color, zorder=6)

    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_aspect("equal")
    if title is not None:
        ax.set_title(title)

    return ax


# ===========================================================================
# Composite building
# ===========================================================================

def build_composite(
    ds: xr.Dataset,
    model_times: np.ndarray,
    comp_info: xr.Dataset,
    scalar_vars: Sequence[str] = ("DBZ", "TH", "QV", "P"),
    vector_vars: Sequence[str] = ("U", "V"),
    derived_vars: Sequence[str] = ("PTE", "UC", "VC"),
    time_coord: str = "time",
    plot_locs: bool = True,
) -> xr.Dataset:
    """
    Build a storm-centred composite from a normalised model Dataset.

    Parameters
    ----------
    ds : xr.Dataset
        Normalised model dataset.
    model_times : np.ndarray
        Model output times (s) to include.
    comp_info : xr.Dataset
        Output of :func:`init_composite`.
    scalar_vars, vector_vars, derived_vars : sequences of str
        Canonical variable names to accumulate.
    time_coord : str
        Name of the time coordinate in *ds*.
    plot_locs : bool
        If ``True``, save a figure of feature-centre locations to disk.

    Returns
    -------
    xr.Dataset of composite-mean arrays on the storm-centred grid.
    """
    ichw = comp_info.attrs["ichw"]
    jchw = comp_info.attrs["jchw"]
    ishw = comp_info.attrs["ishw"]
    jshw = comp_info.attrs["jshw"]
    tracking_varname = comp_info.attrs["tracking_varname"]
    tracking_thresh = comp_info.attrs.get("tracking_thresh", 0.0)
    tracking_extremum = comp_info.attrs["tracking_extremum"]
    tracking_klvl = comp_info.attrs["tracking_klvl"]
    tracking_level = comp_info.attrs["tracking_level"]
    xc1dg = comp_info.coords["xc_g"].values
    yc1dg = comp_info.coords["yc_g"].values
    xckm = comp_info.coords["xc_km"].values
    yckm = comp_info.coords["yc_km"].values
    zeagl1d = comp_info.coords["ze_agl"].values
    zcagl1d = comp_info.coords["zc_agl"].values
    igbgn = comp_info.attrs["igbgn"]
    igend = comp_info.attrs["igend"] - 1
    jgbgn = comp_info.attrs["jgbgn"]
    jgend = comp_info.attrs["jgend"] - 1
    ntimes = model_times.size

    all_scalar_derived = list(dict.fromkeys(list(scalar_vars) + list(derived_vars)))

    varcompdict: dict[str, np.ndarray] = {}
    for name in all_scalar_derived:
        varcompdict[name] = np.zeros((jchw * 2 + 1, ichw * 2 + 1), dtype=np.float64)
    for name in vector_vars:
        varcompdict[name] = np.zeros((jchw * 2 + 2, ichw * 2 + 2), dtype=np.float64)

    if plot_locs:
        figcl, axcl = plt.subplots()
        colour_cycle = [plt.cm.plasma(i) for i in np.linspace(0, 1, ntimes)]

    # Initialise search box at the centre of the sub-domain.
    iref = int((igend - igbgn) / 2.0)
    jref = int((jgend - jgbgn) / 2.0)
    searchboxlims = (
        iref - ishw * 2, iref + ishw * 2 + 1,
        jref - jshw * 2, jref + jshw * 2 + 1,
    )
    guesscoords = (iref, jref)

    nt = 0
    all_var_names = list(
        dict.fromkeys([tracking_varname] + list(scalar_vars) + list(vector_vars))
    )

    for ti, time_s in enumerate(model_times):
        print(f"Model time: {int(time_s):d} s")
        fields = get_fields_at_time(
            ds, float(time_s), all_var_names,
            level=tracking_klvl, time_coord=time_coord,
        )

        trackvar = fields[tracking_varname].values if tracking_varname in fields else None
        if trackvar is None:
            print(
                f"Warning: {tracking_varname!r} not found at "
                f"time {time_s} – skipping."
            )
            continue

        var_extremum, iref, jref = trackfeature(
            trackvar,
            searchboxlims=searchboxlims,
            guesscoords=guesscoords,
            extremum=tracking_extremum,
        )
        guesscoords = (iref, jref)
        searchboxlims = (
            iref - ishw, iref + ishw + 1,
            jref - jshw, jref + jshw + 1,
        )
        print(
            f"{tracking_extremum} {tracking_varname} at z={tracking_level:.1f} m "
            f"= {var_extremum:.2f}  at x,y=({xckm[iref]:.2f}, {yckm[jref]:.2f}) "
            f"i,j=({iref:d},{jref:d})"
        )

        if np.abs(var_extremum) < tracking_thresh:
            continue

        if plot_locs:
            axcl.plot([xckm[iref]], [yckm[jref]], marker="o", color=colour_cycle[ti])

        # Crop scalar and wind fields to composite box.
        js = jref - jchw + jgbgn
        je = jref + jchw + 1 + jgbgn
        is_ = iref - ichw + igbgn
        ie = iref + ichw + 1 + igbgn

        for name in scalar_vars:
            if name in fields:
                varcompdict[name] += fields[name].values[js:je, is_:ie]

        for name in vector_vars:
            if name in fields:
                varcompdict[name] += fields[name].values[js: je + 1, is_: ie + 1]

        if "UC" in derived_vars and "U" in fields:
            u = fields["U"].values
            varcompdict["UC"] += 0.5 * (u[js:je, is_:ie] + u[js:je, is_+1:ie+1])
        if "VC" in derived_vars and "V" in fields:
            v = fields["V"].values
            varcompdict["VC"] += 0.5 * (v[js:je, is_:ie] + v[js+1:je+1, is_:ie])

        if "PTE" in derived_vars and {"P", "TH", "QV"}.issubset(fields):
            try:
                pte = thermo.calpte(
                    fields["P"].values, fields["TH"].values, fields["QV"].values
                )
                varcompdict["PTE"] += pte[js:je, is_:ie]
            except Exception:
                print("Cannot calculate PTE – skipping.")

        nt += 1

    if nt > 1:
        for name in varcompdict:
            varcompdict[name] = varcompdict[name] / nt

    if plot_locs:
        axcl.set_aspect("equal")
        run_label = ds.attrs.get("run_name", "run")
        t0_s, t1_s = int(model_times[0]), int(model_times[-1])
        figcl.savefig(
            f"{run_label}_sfcvortloc_{t0_s:06d}_{t1_s:06d}.png", dpi=200
        )

    xc_comp = np.arange(-ichw, ichw + 1) * comp_info.attrs["dx"] / 1000.0
    yc_comp = np.arange(-jchw, jchw + 1) * comp_info.attrs["dy"] / 1000.0
    xe_comp = np.arange(-ichw - 0.5, ichw + 1.5) * comp_info.attrs["dx"] / 1000.0
    ye_comp = np.arange(-jchw - 0.5, jchw + 1.5) * comp_info.attrs["dy"] / 1000.0
    all_scalar_names = list(dict.fromkeys(list(scalar_vars) + list(derived_vars)))
    data_vars: dict = {}
    for name in all_scalar_names:
        if name in varcompdict:
            data_vars[name] = (["yc_comp", "xc_comp"], varcompdict[name])
    for name in vector_vars:
        if name in varcompdict:
            data_vars[name] = (["ye_comp", "xe_comp"], varcompdict[name])
    return xr.Dataset(
        data_vars,
        coords={
            "xc_comp": xc_comp,
            "yc_comp": yc_comp,
            "xe_comp": xe_comp,
            "ye_comp": ye_comp,
        },
        attrs={"ntimes": nt},
    )


# ===========================================================================
# Observed transect DSD analysis
# ===========================================================================

def _make_logND(ND: np.ndarray) -> np.ma.MaskedArray:
    """Convert ND (m^-3 mm^-1) to masked log10(N / 1000)."""
    logND = np.ma.log10(ND / 1000.0)
    return np.ma.masked_where(logND <= -1.0, logND)


def _plot_obs_transect(
    sample_xlocs: np.ndarray,
    Dl: np.ndarray,
    Dmax_index: int,
    Dmax: float,
    logND: np.ma.MaskedArray,
    D0r: np.ndarray,
    logND_gam: np.ma.MaskedArray | None,
    D0r_gam: np.ndarray | None,
    logND_tmf: np.ma.MaskedArray | None,
    D0r_tmf: np.ndarray | None,
) -> None:
    """Plot observed DSD transect panels (side-effect only)."""
    calc_fits = logND_gam is not None
    nrows = 3 if calc_fits else 1
    fig, axes = plt.subplots(nrows, 1, figsize=(8, 3 * nrows), squeeze=False)
    axes = [a for row in axes for a in row]

    xkm = [sample_xlocs / 1000.0]
    yvals = [Dl[:Dmax_index + 1] * 1000.0]
    _base_params = {
        "type": "pcolor",
        "vlimits": (-1.0, 3.0),
        "clabel": r"log[N ($m^{-3} mm^{-1}$)]",
    }

    panels = [(logND, D0r, not calc_fits)]
    if calc_fits:
        panels += [(logND_gam, D0r_gam, False), (logND_tmf, D0r_tmf, True)]

    for ax, (zfield, d0r, show_cbar) in zip(axes, panels):
        pd = {**_base_params, "plotcbar": show_cbar}
        ax = pm.plotmeteogram(ax, xkm, [zfield.T], [pd], yvals=yvals)
        ax.set_ylim(0.0, Dmax)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(base=2.0))
        ax.set_ylabel("D (mm)")
        ax.plot(xkm[0], d0r * 1000.0, c="k", ls="-", lw=1)

    axes[-1].set_xlabel("x position (km)")


def calc_obs_transect(
    probe_names: Sequence[str],
    probe_nd_list: Sequence[np.ndarray],
    probe_rho_list: Sequence[np.ndarray],
    transect_results: Sequence[xr.Dataset],
    Dr: np.ndarray,
    Dmax: float | None = None,
    calc_fits: bool = True,
    plot_transects: bool = False,
) -> list[xr.Dataset]:
    """
    Compute DSD quantities from observations along probe transects.

    Parameters
    ----------
    probe_names : sequence of str
    probe_nd_list : sequence of np.ndarray
        Per-probe ND arrays, shape ``(ntimes, nbins)``.
    probe_rho_list : sequence of np.ndarray
        Per-probe air-density arrays (kg m^-3), shape ``(ntimes,)``.
    transect_results : sequence of xr.Dataset
        Output of :func:`find_transect_grid_intersections`; the variable
        ``'x_sample'`` provides the sample x-positions.
    Dr : np.ndarray
        Right-edge diameter bins (m).
    Dmax : float or None
        Maximum diameter (m).
    calc_fits : bool
        If ``True``, compute exponential/gamma/TMM DSD fits.
    plot_transects : bool

    Returns
    -------
    list of xr.Dataset, one per probe, with variables ``ND``, ``D0r_obs``,
    and (if *calc_fits*) ``ND_gam``, ``D0r_gam``, ``ND_tmf``, ``D0r_tmf``.
    """
    Dmax, Dmax_index = get_Dmax_index(Dr, Dmax)

    D0r_obs: list[np.ndarray] = []
    ND_obs: list[np.ndarray] = []
    D0r_obs_gam: list[np.ndarray] = []
    D0r_obs_tmf: list[np.ndarray] = []
    ND_obs_gam: list[np.ndarray] = []
    ND_obs_tmf: list[np.ndarray] = []

    for d, _dis_name in enumerate(probe_names):
        sample_xlocs = transect_results[d]["x_sample"].values
        ND = probe_nd_list[d]
        rho = probe_rho_list[d]

        logND_gam = logND_tmf = D0r_gam = D0r_tmf = None
        ND_gamDSD = ND_tmfDSD = None

        if calc_fits:
            _synthbins, exp_DSD, gam_DSD, tmf_DSD, _dis_DSD = dis.calc_DSD(ND.T, rho)

            (ND_gamDSD, _N0_gam, _lamda_gam, _mu_gam, _qr_gam, _Ntr_gam,
             _refl_gam, _D_med_gam, _D_m_gam, _LWC_gam, _rr_gam) = gam_DSD

            (ND_tmfDSD, _N0_tmf, _lamda_tmf, _mu_tmf, _qr_tmf, _Ntr_tmf,
             _refl_tmf, _LWC_tmf, _rr_tmf) = tmf_DSD

            ND_gamDSD = ND_gamDSD[:, :Dmax_index + 1]
            ND_tmfDSD = ND_tmfDSD[:, :Dmax_index + 1]
            logND_gam = _make_logND(ND_gamDSD)
            logND_tmf = _make_logND(ND_tmfDSD)

            D0r_gam = np.array([
                dis.calc_D0_bin(
                    D[:Dmax_index + 1], Dl[:Dmax_index + 1], Dr[:Dmax_index + 1],
                    ND_gamDSD[t, :], bin_width[:Dmax_index + 1],
                )
                for t in range(sample_xlocs.size)
            ])
            D0r_tmf = np.array([
                dis.calc_D0_bin(
                    D[:Dmax_index + 1], Dl[:Dmax_index + 1], Dr[:Dmax_index + 1],
                    ND_tmfDSD[t, :], bin_width[:Dmax_index + 1],
                )
                for t in range(sample_xlocs.size)
            ])
            D0r_obs_gam.append(D0r_gam)
            D0r_obs_tmf.append(D0r_tmf)
            ND_obs_gam.append(ND_gamDSD)
            ND_obs_tmf.append(ND_tmfDSD)

        ND_trunc = ND[:, :Dmax_index + 1]
        logND = _make_logND(ND_trunc)

        D0r = np.array([
            dis.calc_D0_bin(
                D[:Dmax_index + 1], Dl[:Dmax_index + 1], Dr[:Dmax_index + 1],
                ND_trunc[t, :], bin_width[:Dmax_index + 1],
            )
            for t in range(sample_xlocs.size)
        ])
        D0r_obs.append(D0r)
        ND_obs.append(ND_trunc)

        if plot_transects:
            _plot_obs_transect(
                sample_xlocs, Dl, Dmax_index, Dmax,
                logND, D0r, logND_gam, D0r_gam, logND_tmf, D0r_tmf,
            )

    probe_datasets: list[xr.Dataset] = []
    for d in range(len(probe_names)):
        ds_vars: dict = {
            "ND": xr.DataArray(ND_obs[d], dims=["sample_time", "diameter_bin"]),
            "D0r_obs": xr.DataArray(D0r_obs[d], dims=["sample_time"]),
        }
        if calc_fits:
            ds_vars.update({
                "ND_gam": xr.DataArray(ND_obs_gam[d], dims=["sample_time", "diameter_bin"]),
                "D0r_gam": xr.DataArray(D0r_obs_gam[d], dims=["sample_time"]),
                "ND_tmf": xr.DataArray(ND_obs_tmf[d], dims=["sample_time", "diameter_bin"]),
                "D0r_tmf": xr.DataArray(D0r_obs_tmf[d], dims=["sample_time"]),
            })
        probe_datasets.append(xr.Dataset(ds_vars))
    return probe_datasets


# ===========================================================================
# Model transect DSD analysis
# ===========================================================================

def _plot_model_transect(
    sample_xlocs: np.ndarray,
    Dl: np.ndarray,
    Dmax: float,
    Dmax_index: int,
    logNc_bin: np.ma.MaskedArray,
    D0r: np.ndarray,
    logNc_bin_ps: np.ma.MaskedArray | None = None,
    D0r_ps: np.ndarray | None = None,
) -> None:
    """Plot model DSD transect panels (side-effect only)."""
    use_Parsivel_simulator = logNc_bin_ps is not None
    nrows = 2 if use_Parsivel_simulator else 1
    fig, axes = plt.subplots(nrows, 1, figsize=(8, 3 * nrows), squeeze=False)
    axes = [a for row in axes for a in row]

    xkm = [sample_xlocs / 1000.0]
    yvals = [Dl[:Dmax_index + 1] * 1000.0]
    _base_params = {
        "type": "pcolor",
        "vlimits": (-1.0, 3.0),
        "clabel": r"log[N ($m^{-3} mm^{-1}$)]",
    }

    panels = [(logNc_bin, D0r, not use_Parsivel_simulator)]
    if use_Parsivel_simulator:
        panels.append((logNc_bin_ps, D0r_ps, True))

    for ax, (zfield, d0r, show_cbar) in zip(axes, panels):
        pd = {**_base_params, "plotcbar": show_cbar}
        ax = pm.plotmeteogram(ax, xkm, [zfield.T], [pd], yvals=yvals)
        ax.set_ylim(0.0, Dmax)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1.0))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(base=2.0))
        ax.set_ylabel("D (mm)")
        ax.plot(xkm[0], d0r * 1000.0, c="k", ls="-", lw=1)

    axes[-1].set_xlabel("x position (km)")


def interp_model_to_transect(
    probe_names: Sequence[str],
    transect_results: Sequence[xr.Dataset],
    Dr: np.ndarray,
    sampling_interval: float = 60.0,
    sampling_length: float = sampling_length_default,
    sampling_width: float = sampling_width_default,
    add_hail: bool = False,
    use_bins_for_interp: bool = False,
    use_Parsivel_simulator: bool = False,
    Dmax: float | None = None,
    plot_transects: bool = False,
) -> list[xr.Dataset]:
    """
    Interpolate model variables (including bulk DSDs) to simulated
    disdrometer transects and bin to Parsivel diameter bins.

    Parameters
    ----------
    probe_names : sequence of str
    transect_results : sequence of xr.Dataset
        Output of :func:`find_transect_grid_intersections`.  Each Dataset
        must contain variables ``'rhoa'``, ``'qr'``, ``'ntr'``, ``'zr'``,
        ``'alphar'``, etc. at ``all_time``, plus coordinates
        ``'sample_time'`` and ``'all_time'`` and variable
        ``'x_sample'``.
    Dr : np.ndarray
        Right-edge diameter bins (m).
    sampling_interval : float
        Parsivel accumulation interval (s), default 60.
    add_hail : bool
        If ``True``, add hail and graupel contributions to the rain DSD.
        Only valid when *use_bins_for_interp* is ``True``.
    use_bins_for_interp : bool
        If ``True``, discretise the model gamma DSD to Parsivel bins first,
        then interpolate.  If ``False``, use moment-weighted averages.
    use_Parsivel_simulator : bool
        If ``True``, sample the model DSD using the Parsivel simulator.
    Dmax : float or None
        Maximum diameter (m).
    plot_transects : bool

    Returns
    -------
    list of xr.Dataset, one per probe, with variables ``ND``, ``D0r``,
    and (if *use_Parsivel_simulator*) ``ND_ps``, ``D0r_ps``.
    """
    Dmax, Dmax_index = get_Dmax_index(Dr, Dmax)
    sampling_area = sampling_length * sampling_width

    if not use_bins_for_interp and add_hail:
        print(
            "Moment-weighted resampling is not currently compatible with "
            "add_hail=True.  Setting add_hail=False."
        )
        add_hail = False

    D0r_mod: list[np.ndarray] = []
    ND_list: list[np.ndarray] = []
    D0r_mod_ps: list[np.ndarray] = []
    ND_ps_list: list[np.ndarray] = []

    for d, _dis_name in enumerate(probe_names):
        _ds = transect_results[d]
        sampling_times = _ds.coords["sample_time"].values
        all_times = _ds.coords["all_time"].values
        dt = all_times[1:] - all_times[:-1]
        sample_xlocs = _ds["x_sample"].values
        ntimes = _ds.dims["all_time"]

        rhoa = _ds["rhoa"].values
        qr = _ds["qr"].values
        ntr = _ds["ntr"].values
        zr = _ds["zr"].values

        alphar = _ds["alphar"].values
        N0r, _ = dsd.cal_N0(rhoa, qr, ntr, cr, alphar)
        lamdar = dsd.cal_lamda(rhoa, qr, ntr, cr, alphar)

        if add_hail:
            qh = _ds["qh"].values
            nth = _ds["nth"].values
            alphah = _ds["alphah"].values
            rhoh = _ds["rhoh"].values
            ch = rhoh * np.pi / 6.0

            qg = _ds["qg"].values
            ntg = _ds["ntg"].values
            alphag = _ds["alphag"].values
            rhog = _ds["rhog"].values
            cg = rhog * np.pi / 6.0

            N0h, _ = dsd.cal_N0(rhoa, qh, nth, ch, alphah)
            lamdah = dsd.cal_lamda(rhoa, qh, nth, ch, alphah)
            N0g, _ = dsd.cal_N0(rhoa, qg, ntg, cg, alphag)
            lamdag = dsd.cal_lamda(rhoa, qg, ntg, cg, alphag)

        if use_Parsivel_simulator:
            Vtr = np.array([
                dis.assignfallspeed(dis.avg_diameter, rhocorrect=True, rho=rhoa[t])
                for t in range(ntimes)
            ])
            Vtr = Vtr[:, :Dmax_index + 1]

        Nc_bin_tmp = np.empty((np.size(N0r), np.size(D[:Dmax_index + 1])))
        Nc_bin = np.zeros(
            (np.size(np.array(sampling_times)), np.size(D[:Dmax_index + 1]))
        )

        if use_bins_for_interp:
            for index, _ in np.ndenumerate(N0r):
                Nc_bin_tmp[index, :] = (
                    1.0e-3 * N0r[index]
                    * D[:Dmax_index + 1] ** alphar[index]
                    * np.exp(-lamdar[index] * D[:Dmax_index + 1])
                )
                if add_hail:
                    Nc_bin_tmp[index, :] += 1.0e-3 * N0h[index] * (
                        D[:Dmax_index + 1] ** alphah[index]
                        * np.exp(-lamdah[index] * D[:Dmax_index + 1])
                    )
                    Nc_bin_tmp[index, :] += 1.0e-3 * N0g[index] * (
                        D[:Dmax_index + 1] ** alphag[index]
                        * np.exp(-lamdag[index] * D[:Dmax_index + 1])
                    )
        else:
            qr_stimes = np.zeros(np.size(np.array(sampling_times)))
            ntr_stimes = np.zeros_like(qr_stimes)
            zr_stimes = np.zeros_like(qr_stimes)
            N0r_stimes = np.zeros_like(qr_stimes)
            lamdar_stimes = np.zeros_like(qr_stimes)
            alphar_stimes = np.zeros_like(qr_stimes)
            rhoa_stimes = np.zeros_like(qr_stimes)

            # First time: assume constant DSD for the preceding interval.
            Nc_bin_tmp[0, :] = (
                1.0e-3 * N0r[0]
                * D[:Dmax_index + 1] ** alphar[0]
                * np.exp(-lamdar[0] * D[:Dmax_index + 1])
            )

        Nc_bin_tmp = np.ma.masked_invalid(Nc_bin_tmp)

        if use_Parsivel_simulator:
            Nc_bin_tmp_ps = np.empty(
                (np.size(lamdar), np.size(D[:Dmax_index + 1]))
            )
            Nc_bin_ps = np.zeros(
                (np.size(np.array(sampling_times)), np.size(D[:Dmax_index + 1]))
            )
            sample_dict = create_random_gamma_DSD(
                ntr[0], lamdar[0], alphar[0], Vtr[0],
                sampling_length, sampling_width,
                Dl, D, Dr, Dmax=Dmax,
                sampling_interval=sampling_interval,
                remove_margins=True, rhocorrect=True, rho=rhoa[0],
            )
            ND_sample = sample_dict["ND"].values
            pcount_binned_sample = sample_dict["pcount_binned"].values

            if add_hail:
                sample_dict_h = create_random_gamma_DSD(
                    nth[0], lamdah[0], alphah[0], Vtr[0],
                    sampling_length, sampling_width,
                    Dl, D, Dr, Dmax=Dmax,
                    sampling_interval=sampling_interval,
                    remove_margins=True, rhocorrect=True, rho=rhoa[0],
                )
                sample_dict_g = create_random_gamma_DSD(
                    ntg[0], lamdag[0], alphag[0], Vtr[0],
                    sampling_length, sampling_width,
                    Dl, D, Dr, Dmax=Dmax,
                    sampling_interval=sampling_interval,
                    remove_margins=True, rhocorrect=True, rho=rhoa[0],
                )
                ND_sample = (
                    ND_sample
                    + sample_dict_h["ND"].values
                    + sample_dict_g["ND"].values
                )
                pcount_binned_sample = (
                    pcount_binned_sample
                    + sample_dict_h["pcount_binned"].values
                    + sample_dict_g["pcount_binned"].values
                )

            Nc_bin_tmp_ps[0, :] = 1.0e-3 * ND_sample
            Nc_bin_ps[0, :] = Nc_bin_tmp_ps[0, :]

            pcount_binned_samples = []
            for index, _ in np.ndenumerate(lamdar[:-1]):
                sample_dict = create_random_gamma_DSD(
                    ntr[index], lamdar[index], alphar[index], Vtr[index],
                    sampling_length, sampling_width,
                    Dl, D, Dr, Dmax=Dmax,
                    sampling_interval=dt[index],
                    remove_margins=True, rhocorrect=True, rho=rhoa[index],
                )
                ND_sample = sample_dict["ND"].values
                pcount_binned_sample = sample_dict["pcount_binned"].values

                if add_hail:
                    sample_dict_h = create_random_gamma_DSD(
                        nth[index], lamdah[index], alphah[index], Vtr[index],
                        sampling_length, sampling_width,
                        Dl, D, Dr, Dmax=Dmax,
                        sampling_interval=dt[index],
                        remove_margins=True, rhocorrect=True, rho=rhoa[index],
                    )
                    sample_dict_g = create_random_gamma_DSD(
                        ntg[index], lamdag[index], alphag[index], Vtr[index],
                        sampling_length, sampling_width,
                        Dl, D, Dr, Dmax=Dmax,
                        sampling_interval=sampling_interval,
                        remove_margins=True, rhocorrect=True, rho=rhoa[index],
                    )
                    ND_sample = (
                        ND_sample
                        + sample_dict_h["ND"].values
                        + sample_dict_g["ND"].values
                    )
                    pcount_binned_sample = (
                        pcount_binned_sample
                        + sample_dict_h["pcount_binned"].values
                        + sample_dict_g["pcount_binned"].values
                    )

                pcount_binned_samples.append(sample_dict["pcount_binned"].values)
                Nc_bin_tmp_ps[index, :] = 1.0e-3 * ND_sample

            pcount_binned_samples = np.array(pcount_binned_samples)
            Nc_bin_tmp_ps = np.ma.masked_invalid(Nc_bin_tmp_ps)

        else:
            # First time: assume constant DSD for the preceding interval.
            sample_indices = np.searchsorted(all_times, sampling_times, side="left")
            Nc_bin[0, :] = Nc_bin_tmp[sample_indices[0], :]

        sample_indices = np.searchsorted(all_times, sampling_times, side="left")

        for s, sample_index in enumerate(sample_indices[:-1]):
            sample_index_end = sample_indices[s + 1]
            current_slice = slice(sample_index, sample_index_end, None)

            if use_bins_for_interp:
                Nc_bin[s + 1, :] = (
                    np.sum(
                        Nc_bin_tmp[current_slice, :] * dt[current_slice, None],
                        axis=0,
                    )
                    / sampling_interval
                )
            else:
                qr_stimes[s + 1] = (
                    np.sum(qr[current_slice] * dt[current_slice]) / sampling_interval
                )
                ntr_stimes[s + 1] = (
                    np.sum(ntr[current_slice] * dt[current_slice]) / sampling_interval
                )
                zr_stimes[s + 1] = (
                    np.sum(zr[current_slice] * dt[current_slice]) / sampling_interval
                )
                rhoa_stimes[s + 1] = (
                    np.sum(rhoa[current_slice] * dt[current_slice]) / sampling_interval
                )
                alphar_stimes[s + 1] = dualpol.solve_alpha_iter(
                    rhoa_stimes[s + 1], mur,
                    qr_stimes[s + 1], ntr_stimes[s + 1],
                    zr_stimes[s + 1], rhorcst,
                )
                N0r_stimes[s + 1], _ = dsd.cal_N0(
                    rhoa_stimes[s + 1], qr_stimes[s + 1],
                    ntr_stimes[s + 1], cr, alphar_stimes[s + 1],
                )
                lamdar_stimes[s + 1] = dsd.cal_lamda(
                    rhoa_stimes[s + 1], qr_stimes[s + 1],
                    ntr_stimes[s + 1], alphar_stimes[s + 1],
                )
                Nc_bin[s + 1, :] = (
                    1.0e-3
                    * N0r_stimes[s + 1]
                    * D[:Dmax_index + 1] ** alphar_stimes[s + 1]
                    * np.exp(-lamdar_stimes[s + 1] * D[:Dmax_index + 1])
                )

            if use_Parsivel_simulator:
                Vtr_mean = (
                    np.sum(
                        Vtr[current_slice, :] * dt[current_slice, None], axis=0
                    )
                    / sampling_interval
                )
                _samp_vols_D = calc_sampling_volumes_D(
                    Vtr_mean, Dr, Dmax, sampling_interval, sampling_area
                )
                pcount_binned = np.sum(
                    pcount_binned_samples[current_slice], axis=0
                )
                Nc_bin_ps[s + 1, :] = 1.0e-3 * calc_ND(
                    pcount_binned, _samp_vols_D, Dr, Dl, Dmax
                )

        logNc_bin = np.ma.masked_where(
            np.log10(Nc_bin) <= -1.0, np.log10(Nc_bin)
        )

        logNc_bin_ps = None
        if use_Parsivel_simulator:
            Nc_bin_ps = np.ma.masked_invalid(Nc_bin_ps)
            logNc_bin_ps = np.ma.masked_where(
                np.log10(Nc_bin_ps) <= -1.0, np.log10(Nc_bin_ps)
            )

        D0r = np.array([
            dis.calc_D0_bin(
                D[:Dmax_index + 1], Dl[:Dmax_index + 1], Dr[:Dmax_index + 1],
                Nc_bin[t, :], bin_width[:Dmax_index + 1],
            )
            for t in range(sampling_times.size)
        ])
        D0r_mod.append(D0r)
        ND_list.append(Nc_bin)

        if use_Parsivel_simulator:
            D0r_ps = np.array([
                dis.calc_D0_bin(
                    D[:Dmax_index + 1], Dl[:Dmax_index + 1], Dr[:Dmax_index + 1],
                    Nc_bin_ps[t, :], bin_width[:Dmax_index + 1],
                )
                for t in range(sampling_times.size)
            ])
            D0r_mod_ps.append(D0r_ps)
            ND_ps_list.append(Nc_bin_ps)

        if plot_transects:
            _plot_model_transect(
                sample_xlocs, Dl, Dmax, Dmax_index,
                logNc_bin, D0r,
                logNc_bin_ps if use_Parsivel_simulator else None,
                D0r_ps if use_Parsivel_simulator else None,
            )

    probe_datasets: list[xr.Dataset] = []
    for d in range(len(probe_names)):
        ds_vars: dict = {
            "ND": xr.DataArray(ND_list[d], dims=["sample_time", "diameter_bin"]),
            "D0r": xr.DataArray(D0r_mod[d], dims=["sample_time"]),
        }
        if use_Parsivel_simulator:
            ds_vars.update({
                "ND_ps": xr.DataArray(ND_ps_list[d], dims=["sample_time", "diameter_bin"]),
                "D0r_ps": xr.DataArray(D0r_mod_ps[d], dims=["sample_time"]),
            })
        probe_datasets.append(xr.Dataset(ds_vars))
    return probe_datasets


# ===========================================================================
# ARPS ensemble helpers
# ===========================================================================

def get_ARPS_member_dir_and_prefix(member: int, cycle: str) -> tuple[str, str]:
    """
    Return ``(member_dir, member_prefix)`` for a given ensemble member.

    Parameters
    ----------
    member : int
        Member number; 0 is interpreted as the ensemble mean.
    cycle : str
        ``'posterior'`` or ``'prior'``.

    Returns
    -------
    (member_dir, member_prefix)
    """
    if member == 0:
        if "posterior" in cycle:
            return "ENamean", "enmean"
        return "ENfmean", "efmean"
    if "posterior" in cycle:
        return f"EN{member:03d}", f"ena{member:03d}"
    return f"ENF{member:03d}", f"enf{member:03d}"


def read_ARPS_ensemble(
    f,
    member_list: list[int],
    f_args: tuple = (),
    f_kwargs: dict | None = None,
    iterate_over: str = "member_list",
    process_parallel: bool = True,
    n_jobs: int = 5,
    verbose: int = 0,
) -> list:
    """
    Read multiple ARPS ensemble members in serial or parallel (joblib).

    Parameters
    ----------
    f : callable
        Function to call for each member.  It must have a positional argument
        whose name matches *iterate_over*.
    member_list : list of int
        Member numbers to process.
    f_args : tuple
        Positional arguments for *f*.
    f_kwargs : dict or None
        Keyword arguments for *f*.
    iterate_over : str
        Name of the argument in *f* to iterate over.
    process_parallel : bool
    n_jobs : int
        Number of parallel jobs.
    verbose : int
        joblib verbosity level.

    Returns
    -------
    list of results from *f*.
    """
    if f_kwargs is None:
        f_kwargs = {}

    f_argnames = inspect.getfullargspec(f).args
    try:
        arg_to_iterate = f_argnames.index(iterate_over)
    except ValueError:
        print(f"{iterate_over!r} not found in the function argument list – stopping.")
        return []

    def _g(val):
        new_args = tuple(
            val if i == arg_to_iterate else f_args[i]
            for i in range(len(f_args))
        )
        return f(*new_args, **f_kwargs)

    if process_parallel:
        return Parallel(n_jobs=n_jobs, verbose=verbose)(
            delayed(_g)(val) for val in member_list
        )
    return [_g(val) for val in member_list]


def read_ARPS_member_data(
    basedir: str,
    expname: str,
    member: int,
    cycle: str,
    fileformat: str,
    time_range: Sequence,
    tintv_mean: float,
    varnames: Sequence[str],
    filetype: str = "history",
    ibgn: int | None = None,
    iend: int | None = None,
    jbgn: int | None = None,
    jend: int | None = None,
    klvls: list[int] | None = None,
    nproc_x: int = 1,
    nproc_y: int = 1,
    dump_nc_output: bool = True,
    ncdir: str | None = None,
    fileprefix: str = "",
    ds_grid: xr.Dataset | None = None,
    datetime_range: Sequence | None = None,
    x: np.ndarray | None = None,
    y: np.ndarray | None = None,
    mid_diameters: np.ndarray | None = None,
) -> xr.Dataset:
    """
    Read data for a single ARPS ensemble member across multiple times.

    Parameters
    ----------
    basedir : str
        Base directory of the ensemble.
    expname : str
        Experiment name (sub-directory of *basedir*).
    member : int
    cycle : str
        ``'posterior'`` or ``'prior'``.
    fileformat : str
        ARPS file format string.
    time_range : sequence
        Model times to read.
    tintv_mean : float
        Interval between analysis (mean) times.
    varnames : sequence of str
        Variable names to read.
    filetype : str
        ARPS file type, default ``'history'``.
    ibgn, iend, jbgn, jend : int or None
        Optional spatial sub-domain indices.
    klvls : list of int or None
        Vertical levels to read.
    nproc_x, nproc_y : int
        Decomposition dimensions.
    dump_nc_output : bool
        If ``True``, write the read data to a netCDF file via
        :func:`dump_ARPS_xyslice_nc`.
    ncdir : str or None
        Output directory for netCDF files.
    fileprefix : str
        Prefix prepended to output filenames.
    ds_grid : xr.Dataset or None
        Grid dataset providing coordinate arrays and metadata attributes.
    datetime_range : sequence or None
        Datetime objects matching *time_range* for the netCDF time axis.
    x, y : np.ndarray or None
        1-D coordinate arrays used when *ds_grid* is not provided.
    mid_diameters : np.ndarray or None
        Diameter bin centres for DSD output.

    Returns
    -------
    xr.Dataset with all read variables on the model grid, or an empty
    Dataset if coordinate information was not available.
    """
    if not _HAS_PYCRMTOOLS:
        raise ImportError(
            "pyCRMtools is required for read_ARPS_member_data(). "
            "Install it or use xarray-based loading instead."
        )
    if klvls is None:
        klvls = [1]

    print(f"Loading member #{member:d}")
    vardict_list = []

    for time in time_range:
        print("Loading time", time)
        cycle_temp = "posterior" if time % tintv_mean != 0 else "prior"
        member_dir, member_prefix = get_ARPS_member_dir_and_prefix(member, cycle_temp)
        member_absdir = os.path.join(basedir, expname, member_dir)
        filepath = arps_read.get_file_path(
            member_absdir, member_prefix, fileformat, time=time, filetype="history"
        )
        print(filepath)
        vardict = arps_read.read_hdfvars(
            filepath, varnames,
            ibgn=ibgn, jbgn=jbgn, iend=iend, jend=jend,
            klvls=klvls, nproc_x=nproc_x, nproc_y=nproc_y,
        )
        vardict_list.append(vardict)

    # Build xr.Dataset from the list of vardicts when coordinate info is available.
    has_coords = (ds_grid is not None) or (x is not None and y is not None)
    has_indices = all(v is not None for v in (ibgn, iend, jbgn, jend))
    if has_coords and has_indices:
        if ds_grid is not None:
            xc_arr = np.asarray(ds_grid.coords["xc"].values).ravel()
            yc_arr = np.asarray(ds_grid.coords["yc"].values).ravel()
            xe_arr = np.asarray(ds_grid.coords["xe"].values).ravel()
            ye_arr = np.asarray(ds_grid.coords["ye"].values).ravel()
        else:
            xc_arr = yc_arr = x
            xe_arr = ye_arr = y

        coord_dict: dict = {
            "time": datetime_range,
            "yc": ("yc", yc_arr),
            "xc": ("xc", xc_arr),
            "ye": ("ye", ye_arr),
            "xe": ("xe", xe_arr),
        }
        vardict_combined = CRMutils.make_dict_of_lists(vardict_list)
        for varname, var in vardict_combined.items():
            var_arr = np.array(var).T.squeeze()
            var_arr = np.rollaxis(var_arr, 2, 0)
            if varname == "u":
                var_arr_patch = var_arr[:, jbgn: jend + 1, ibgn: iend + 2]
                vardict_combined[varname] = (["time", "yc", "xe"], var_arr_patch)
            elif varname == "v":
                var_arr_patch = var_arr[:, jbgn: jend + 2, ibgn: iend + 1]
                vardict_combined[varname] = (["time", "ye", "xc"], var_arr_patch)
            else:
                var_arr_patch = var_arr[:, jbgn: jend + 1, ibgn: iend + 1]
                vardict_combined[varname] = (["time", "yc", "xc"], var_arr_patch)

        var_ds = xr.Dataset(vardict_combined, coords=coord_dict)
        if ds_grid is not None:
            for attr in (
                "nx_full", "ny_full", "dx", "dy",
                "ctrlat", "ctrlon", "trulat1", "trulat2", "trulon",
            ):
                if attr in ds_grid.attrs:
                    var_ds.attrs[attr] = ds_grid.attrs[attr]
    else:
        var_ds = xr.Dataset()

    if dump_nc_output and var_ds.data_vars:
        output_prefix = fileprefix + member_prefix
        dump_ARPS_xyslice_nc(
            var_ds=var_ds,
            ncdir=ncdir,
            member_prefix=output_prefix,
            mid_diameters=mid_diameters,
        )

    return var_ds


def dump_ARPS_xyslice_nc(
    var_ds: xr.Dataset,
    ncdir: str | None = None,
    member_prefix: str | None = None,
    mid_diameters: np.ndarray | None = None,
) -> None:
    """
    Add derived DSD fields to an ARPS xy-slice Dataset and write to netCDF.

    Parameters
    ----------
    var_ds : xr.Dataset
        Dataset produced by :func:`read_ARPS_member_data`, already containing
        all raw variables on the model grid with proper coordinates.
    ncdir : str or None
        Output directory.
    member_prefix : str or None
        File name prefix.
    mid_diameters : np.ndarray or None
        Diameter bin centres; when provided, the full binned DSD is computed
        and added as ``ND`` and ``logND``.
    """
    if not _HAS_PYCRMTOOLS:
        raise ImportError(
            "pyCRMtools is required for dump_ARPS_xyslice_nc()."
        )

    rhor = 1000.0
    _cr = np.pi / 6.0 * rhor
    var_ds["rho"] = thermo.calrho(var_ds["p"], var_ds["pt"], var_ds["qv"])
    var_ds["alphar"] = dsd.solve_alpha(
        var_ds["rho"], _cr, var_ds["qr"], var_ds["nr"], var_ds["zr"]
    )
    var_ds["N0r"] = dsd.calc_N0_gamma(
        var_ds["rho"], var_ds["qr"], var_ds["nr"], _cr, var_ds["alphar"]
    )
    var_ds["lamdar"] = dsd.calc_lamda_gamma(
        var_ds["rho"], var_ds["qr"], var_ds["nr"], _cr, var_ds["alphar"]
    )

    if mid_diameters is not None:
        _mid_d, N0r_da, lamdar_da, alphar_da = xr.broadcast(
            mid_diameters, var_ds["N0r"], var_ds["lamdar"], var_ds["alphar"]
        )
        _mid_d = _mid_d.transpose("time", "diameter_bin", "yc", "xc")
        N0r_da = N0r_da.transpose("time", "diameter_bin", "yc", "xc")
        lamdar_da = lamdar_da.transpose("time", "diameter_bin", "yc", "xc")
        alphar_da = alphar_da.transpose("time", "diameter_bin", "yc", "xc")

        ND_model = dsd.calc_binned_DSD_from_params(N0r_da, lamdar_da, alphar_da, _mid_d)
        ND_model = ND_model.fillna(0.0)
        logND_model = np.log10(ND_model)
        logND_model = logND_model.where(logND_model > -np.inf)

        var_ds["ND"] = ND_model
        var_ds["logND"] = logND_model

    filename = f"{member_prefix}_fields.nc"
    filepath = os.path.join(ncdir, filename)
    var_ds.to_netcdf(filepath)


# ===========================================================================
# Parsivel-simulator transect sampling from raw model dataset
# ===========================================================================

def sample_model_PSD_along_transect(
    transect_ds: xr.Dataset,
    ds: xr.Dataset,
    grid_idx_ds: xr.Dataset,
    Dmax: float | None = None,
    level: int = 0,
    sampling_length: float = sampling_length_default,
    sampling_width: float = sampling_width_default,
    verbose: bool = False,
) -> xr.Dataset:
    """
    Sample the model PSD along a transect using the Parsivel simulator.

    Single-probe, Parsivel-simulator-only replacement for
    :func:`interp_model_to_transect`.  Model fields are extracted directly
    from *ds* at each transect grid-cell rather than being read from a
    pre-built transect Dataset.  The shape parameter *alphar* is derived
    from the ``zr`` reflectivity field via :func:`dualpara.solve_alpha_iter`.

    Parameters
    ----------
    transect_ds : xr.Dataset
        Output of :func:`find_transect_grid_intersections`.
        Must have coords ``all_time`` and ``sample_time``.
    ds : xr.Dataset
        Normalised model Dataset with canonical field names
        (``rhoa``, ``qr``, ``ntr``, ``zr``).
        Scalar fields are expected to have dimensions
        ``(time, [vertical,] yc, xc)``.
    grid_idx_ds : xr.Dataset
        Output of :func:`get_transect_grid_indices`.
        Variables: ``i_idx``, ``j_idx``, optionally ``time_idx``.
        Coord: ``all_time``.
    Dmax : float or None
        Maximum diameter (m).  Uses the full bin range when ``None``.
    level : int
        Vertical level index (0 = lowest model level).

    Returns
    -------
    xr.Dataset
        Variables:

        ``Nc_bin_ps`` – (sample_time, diameter_bin) number concentration

        Coord: ``sample_time``.
    """
    sampling_area = sampling_length * sampling_width
    Dmax, Dmax_index = get_Dmax_index(Dr, Dmax)

    sampling_times = transect_ds.coords["sample_time"].values
    all_times = transect_ds.coords["all_time"].values

    # Time deltas between consecutive all_time points (seconds).
    dt = all_times[1:] - all_times[:-1]
    # Per-sample-interval durations (seconds), implicit in sampling_times.
    sampling_dt = sampling_times[1:] - sampling_times[:-1]

    i_idx_arr = grid_idx_ds["i_idx"].values.astype(int)
    j_idx_arr = grid_idx_ds["j_idx"].values.astype(int)
    has_time_idx = "time_idx" in grid_idx_ds
    if has_time_idx:
        time_idx_arr = grid_idx_ds["time_idx"].values.astype(int)

    ntimes = len(all_times)
    time_coord = "time"

    # ------------------------------------------------------------------
    # Extract model scalars for every all_time point using bulk numpy
    # indexing to avoid per-step xarray overhead.
    # ------------------------------------------------------------------

    def _get_field_numpy(name: str) -> np.ndarray:
        """Return a 3-D (time, y, x) or 2-D (y, x) numpy array for *name*."""
        vals = ds[name].values
        if vals.ndim == 4:          # (time, level, y, x)
            return vals[:, level, :, :]
        if vals.ndim == 3 and time_coord in ds[name].dims:
            return vals             # (time, y, x)
        if vals.ndim == 3:          # (level, y, x) — no time dimension
            return vals[level]      # → (y, x)
        return vals                 # (y, x)

    def _index_field(arr: np.ndarray, tidx: np.ndarray | int) -> np.ndarray:
        """Advanced-index *arr* to extract the ntimes-long 1-D transect."""
        if arr.ndim == 3:
            return arr[tidx, j_idx_arr, i_idx_arr]
        # 2-D (y, x) — same spatial location for all times
        return arr[j_idx_arr, i_idx_arr]

    if has_time_idx:
        tidx = time_idx_arr
    elif time_coord in ds.dims:
        tidx = np.zeros(ntimes, dtype=int)
    else:
        tidx = None  # no time dimension

    # TODO: make this more general by leveraging the variable mapping in model_config.py
    # Right now, this works for CM1 model output.
    # EDIT: changed this to use the canonical field names in the model dataset, so it should work
    # assuming that the model dataset has been "normalized" to have the canonical field names
    # (rhoa, qr, ntr, zr).
    if "rhoa" in ds:
        rhoa_arr = _index_field(_get_field_numpy("rhoa"), tidx)
    elif "rho" in ds:
        rhoa_arr = _index_field(_get_field_numpy("rho"), tidx)
    else:
        # Derive air density from pressure, potential temperature, and water vapour.
        p_arr  = _index_field(_get_field_numpy("p"), tidx)
        pt_arr = _index_field(_get_field_numpy("th"),  tidx)
        qv_arr = _index_field(_get_field_numpy("qv"),  tidx)
        rhoa_arr = np.array([
            float(thermo.calrho(p_arr[n], pt_arr[n], qv_arr[n]))
            for n in range(ntimes)
        ])

    qr_arr  = _index_field(_get_field_numpy("qr"),  tidx)
    ntr_arr = _index_field(_get_field_numpy("ntr"), tidx)
    zr_arr  = _index_field(_get_field_numpy("zr"), tidx)

    # ------------------------------------------------------------------
    # Derive alphar from reflectivity using the iterative solver.
    # ------------------------------------------------------------------
    alphar_arr = np.array([
        dualpol.solve_alpha_iter(
            rhoa_arr[n], mur, qr_arr[n], ntr_arr[n], zr_arr[n], rhorcst
        ).squeeze()
        for n in range(ntimes)
    ])
    # print("Shape of alphar_arr:", np.shape(alphar_arr))

    # ------------------------------------------------------------------
    # Derived DSD parameters.
    # ------------------------------------------------------------------
    N0r = dsd.calc_N0_gamma(rhoa_arr, qr_arr, ntr_arr, cr, alphar_arr)
    lamdar = dsd.calc_lamda_gamma(rhoa_arr, qr_arr, ntr_arr, cr, alphar_arr)
    # print("Shape of rhoa, qr, ntr, zr:", np.shape(rhoa_arr), np.shape(qr_arr), np.shape(ntr_arr), np.shape(zr_arr))
    # print("Shape of lamdar:", np.shape(lamdar))

    # ------------------------------------------------------------------
    # Terminal-fall speeds for each all_time point.
    # calc_empirical_fallspeed with a 1-D rho Series returns (ntimes, nbins).
    # ------------------------------------------------------------------
    Vtr = pips.calc_empirical_fallspeed(
        D * 1000., correct_rho=True, rho=rhoa_arr
    )
    Vtr = Vtr[:, :Dmax_index + 1]

    # ------------------------------------------------------------------
    # Parsivel simulator: one call per sub-interval between all_times.
    # ------------------------------------------------------------------
    pcount_binned_samples = []
    for n in range(ntimes - 1):
        # print(f"lamdar[{n}]: {lamdar[n]}")
        sample_dict = create_random_gamma_DSD(
            ntr_arr[n], lamdar[n], alphar_arr[n], Vtr[n],
            sampling_length, sampling_width,
            Dl, D, Dr, Dmax=Dmax,
            sampling_interval=float(dt[n]),
            remove_margins=True, rhocorrect=True, rho=rhoa_arr[n],
            verbose=verbose
        )
        if sample_dict is not None:
            pcount_binned_samples.append(sample_dict["pcount_binned"].values)
        else:
            pcount_binned_samples.append(np.zeros(Dmax_index + 1))

    pcount_binned_samples = np.array(pcount_binned_samples)  # (ntimes-1, nbins)

    # ------------------------------------------------------------------
    # Combine sub-interval counts into per-sample-interval concentrations.
    # ------------------------------------------------------------------
    nsamples = len(sampling_times)
    Nc_bin_ps = np.zeros((nsamples, Dmax_index + 1))

    sample_indices = np.searchsorted(all_times, sampling_times, side="left")

    for s, sample_index in enumerate(sample_indices[:-1]):
        sample_index_end = sample_indices[s + 1]
        current_slice = slice(sample_index, sample_index_end, None)

        Vtr_mean = (
            np.sum(Vtr[current_slice, :] * dt[current_slice, None], axis=0)
            / sampling_dt[s]
        )
        _samp_vols_D = calc_sampling_volumes_D(
            Vtr_mean, Dr, Dmax, sampling_dt[s], sampling_area
        )
        pcount_binned_total = np.sum(pcount_binned_samples[current_slice], axis=0)
        Nc_bin_ps[s + 1, :] = 1.0e-3 * calc_ND(
            pcount_binned_total, _samp_vols_D, Dr, Dl, Dmax
        )

    Nc_bin_ps = np.ma.masked_invalid(Nc_bin_ps)

    return xr.Dataset(
        {
            "Nc_bin_ps": xr.DataArray(
                Nc_bin_ps, dims=["sample_time", "diameter_bin"]
            ),
        },
        coords={"sample_time": sampling_times},
    )



## Modified version of sample_model_PSD_along_transect function, but for hail
def sample_model_PSD_along_transect_hail(
    transect_ds: xr.Dataset,
    ds: xr.Dataset,
    grid_idx_ds: xr.Dataset,
    Dmax: float | None = None,
    level: int = 0,
    sampling_length: float = sampling_length_default,
    sampling_width: float = sampling_width_default,
    verbose: bool = False,
) -> xr.Dataset:
    """
    Sample the model PSD along a transect using the Parsivel simulator.

    Single-probe, Parsivel-simulator-only replacement for
    :func:`interp_model_to_transect`.  Model fields are extracted directly
    from *ds* at each transect grid-cell rather than being read from a
    pre-built transect Dataset.  The shape parameter *alphar* is derived
    from the ``zr`` reflectivity field via :func:`dualpara.solve_alpha_iter`.

    Parameters
    ----------
    transect_ds : xr.Dataset
        Output of :func:`find_transect_grid_intersections`.
        Must have coords ``all_time`` and ``sample_time``.
    ds : xr.Dataset
        Normalised model Dataset with canonical field names
        (``rhoa``, ``qr``, ``ntr``, ``zr``).
        Scalar fields are expected to have dimensions
        ``(time, [vertical,] yc, xc)``.
    grid_idx_ds : xr.Dataset
        Output of :func:`get_transect_grid_indices`.
        Variables: ``i_idx``, ``j_idx``, optionally ``time_idx``.
        Coord: ``all_time``.
    Dmax : float or None
        Maximum diameter (m).  Uses the full bin range when ``None``.
    level : int
        Vertical level index (0 = lowest model level).

    Returns
    -------
    xr.Dataset
        Variables:

        ``Nc_bin_ps`` – (sample_time, diameter_bin) number concentration

        Coord: ``sample_time``.
    """
    Dmax, Dmax_index = get_Dmax_index(Dr, Dmax)
    sampling_area = sampling_length * sampling_width
    # print("Dmax, Dmax_index, sampling_area:", Dmax, Dmax_index, sampling_area)

    sampling_times = transect_ds.coords["sample_time"].values
    all_times = transect_ds.coords["all_time"].values

    # Time deltas between consecutive all_time points (seconds).
    dt = all_times[1:] - all_times[:-1]
    # print("dt:", dt)
    # Per-sample-interval durations (seconds), implicit in sampling_times.
    sampling_dt = sampling_times[1:] - sampling_times[:-1]
    # print("sampling_dt:", sampling_dt)

    i_idx_arr = grid_idx_ds["i_idx"].values.astype(int)
    j_idx_arr = grid_idx_ds["j_idx"].values.astype(int)
    has_time_idx = "time_idx" in grid_idx_ds
    if has_time_idx:
        time_idx_arr = grid_idx_ds["time_idx"].values.astype(int)

    ntimes = len(all_times)
    time_coord = "time"

    # ------------------------------------------------------------------
    # Extract model scalars for every all_time point using bulk numpy
    # indexing to avoid per-step xarray overhead.
    # ------------------------------------------------------------------

    def _get_field_numpy(name: str) -> np.ndarray:
        """Return a 3-D (time, y, x) or 2-D (y, x) numpy array for *name*."""
        vals = ds[name].values
        if vals.ndim == 4:          # (time, level, y, x)
            return vals[:, level, :, :]
        if vals.ndim == 3 and time_coord in ds[name].dims:
            return vals             # (time, y, x)
        if vals.ndim == 3:          # (level, y, x) — no time dimension
            return vals[level]      # → (y, x)
        return vals                 # (y, x)

    def _index_field(arr: np.ndarray, tidx: np.ndarray | int) -> np.ndarray:
        """Advanced-index *arr* to extract the ntimes-long 1-D transect."""
        if arr.ndim == 3:
            return arr[tidx, j_idx_arr, i_idx_arr]
        # 2-D (y, x) — same spatial location for all times
        return arr[j_idx_arr, i_idx_arr]

    if has_time_idx:
        tidx = time_idx_arr
    elif time_coord in ds.dims:
        tidx = np.zeros(ntimes, dtype=int)
    else:
        tidx = None  # no time dimension

    # TODO: make this more general by leveraging the variable mapping in model_config.py
    # Right now, this works for CM1 model output.
    # EDIT: changed this to use the canonical field names in the model dataset, so it should work
    # assuming that the model dataset has been "normalized" to have the canonical field names
    # (rhoa, qr, ntr, zr).
    if "rhoa" in ds:
        rhoa_arr = _index_field(_get_field_numpy("rhoa"), tidx)
    elif "rho" in ds:
        rhoa_arr = _index_field(_get_field_numpy("rho"), tidx)
    else:
        # Derive air density from pressure, potential temperature, and water vapour.
        p_arr  = _index_field(_get_field_numpy("p"), tidx)
        pt_arr = _index_field(_get_field_numpy("th"),  tidx)
        qv_arr = _index_field(_get_field_numpy("qv"),  tidx)
        rhoa_arr = np.array([
            float(thermo.calrho(p_arr[n], pt_arr[n], qv_arr[n]))
            for n in range(ntimes)
        ])

    #qr_arr  = _index_field(_get_field_numpy("qr"),  tidx)
    #ntr_arr = _index_field(_get_field_numpy("crw"), tidx)
    #zr_arr  = _index_field(_get_field_numpy("zrw"), tidx)

    qhl_arr = _index_field(_get_field_numpy("qh"), tidx)
    nthl_arr = _index_field(_get_field_numpy("nth"), tidx)
    zhl_arr  = _index_field(_get_field_numpy("zh"), tidx)

    vhl_arr = _index_field(_get_field_numpy("vh"), tidx)

    # calculate density of hail
    rhohl_arr = qhl_arr / vhl_arr
    # For safety, set any bad values of rhohl_arr to 900 kg/m^3
    rhohl_arr = np.where(~np.isfinite(rhohl_arr), 900, rhohl_arr)
    chl_arr = np.pi / 6. * rhohl_arr

    # ------------------------------------------------------------------
    # Derive alphahl from reflectivity using the iterative solver.
    # ------------------------------------------------------------------
    alphahl_arr = np.array([
        dualpol.solve_alpha_iter(
            rhoa_arr[n], muhl, qhl_arr[n], nthl_arr[n], zhl_arr[n], rhohl_arr[n]
        ).squeeze()
        for n in range(ntimes)
    ])
    # print("Shape of alphar_arr:", np.shape(alphar_arr))

    # ------------------------------------------------------------------
    # Derived DSD parameters.
    # ------------------------------------------------------------------
    N0hl = dsd.calc_N0_gamma(rhoa_arr, qhl_arr, nthl_arr, chl_arr, alphahl_arr)
    lamdahl = dsd.calc_lamda_gamma(rhoa_arr, qhl_arr, nthl_arr, chl_arr, alphahl_arr)
    # print("Shape of rhoa, qr, ntr, zr:", np.shape(rhoa_arr), np.shape(qr_arr), np.shape(ntr_arr), np.shape(zr_arr))
    # print("Shape of lamdar:", np.shape(lamdar))

    # ------------------------------------------------------------------
    # Terminal-fall speeds for each all_time point.
    # calc_empirical_fallspeed_hail with a 1-D rho Series returns (ntimes, nbins).
    # ------------------------------------------------------------------
    Vthl = pips.calc_empirical_fallspeed_hail(
        D, rho_h=rhohl_arr, correct_rho=True, rho=rhoa_arr
    )
    Vthl = Vthl[:, :Dmax_index + 1]

    # ------------------------------------------------------------------
    # Parsivel simulator: one call per sub-interval between all_times.
    # ------------------------------------------------------------------
    pcount_binned_samples = []
    for n in range(ntimes - 1):
        # print(f"lamdar[{n}]: {lamdar[n]}")
        sample_dict = create_random_gamma_DSD(
            nthl_arr[n], lamdahl[n], alphahl_arr[n], Vthl[n],
            sampling_length, sampling_width,
            Dl, D, Dr, Dmax=Dmax,
            sampling_interval=float(dt[n]),
            remove_margins=True, rhocorrect=True, rho=rhoa_arr[n], rhohl=rhohl_arr[n],
            verbose=verbose
        )
        if sample_dict is not None:
            pcount_binned_samples.append(sample_dict["pcount_binned"].values)
            # print("sample_dict['pcount_binned'].values:", sample_dict["pcount_binned"].values)
        else:
            pcount_binned_samples.append(np.zeros(Dmax_index + 1))

    pcount_binned_samples = np.array(pcount_binned_samples)  # (ntimes-1, nbins)

    # ------------------------------------------------------------------
    # Combine sub-interval counts into per-sample-interval concentrations.
    # ------------------------------------------------------------------
    nsamples = len(sampling_times)
    Nc_bin_ps = np.zeros((nsamples, Dmax_index + 1))

    sample_indices = np.searchsorted(all_times, sampling_times, side="left")

    for s, sample_index in enumerate(sample_indices[:-1]):
        sample_index_end = sample_indices[s + 1]
        current_slice = slice(sample_index, sample_index_end, None)

        Vthl_mean = (
            np.sum(Vthl[current_slice, :] * dt[current_slice, None], axis=0)
            / sampling_dt[s]
        )
        _samp_vols_D = calc_sampling_volumes_D(
            Vthl_mean, Dr, Dmax, sampling_dt[s], sampling_area
        )
        # print(f"_samp_vols_D: {_samp_vols_D}")
        pcount_binned_total = np.sum(pcount_binned_samples[current_slice], axis=0)
        # print(f"pcount_binned_total: {pcount_binned_total}")
        Nc_bin_ps[s + 1, :] = 1.0e-3 * calc_ND(
            pcount_binned_total, _samp_vols_D, Dr, Dl, Dmax
        )
        # print(f"Nc_bin_ps[s + 1, :]: {Nc_bin_ps[s + 1, :]}")
        # print(f"np.log10(Nc_bin_ps[s + 1, :]): {np.log10(Nc_bin_ps[s + 1, :])}")

    Nc_bin_ps = np.ma.masked_invalid(Nc_bin_ps)

    return xr.Dataset(
        {
            "Nc_bin_ps": xr.DataArray(
                Nc_bin_ps, dims=["sample_time", "diameter_bin"]
            ),
        },
        coords={"sample_time": sampling_times},
    )



def calc_model_PSD_along_transect(
    transect_ds: xr.Dataset,
    ds: xr.Dataset,
    grid_idx_ds: xr.Dataset,
    Dmax: float | None = None,
    level: int = 0,
) -> xr.Dataset:
    """
    Compute the raw model gamma PSD along a transect, binned into Parsivel bins.

    Unlike :func:`sample_model_PSD_along_transect`, this function does not run
    the Parsivel simulator.  Instead it uses the predicted model moments
    (``qr``, ``ntr``, ``zr``) to derive the gamma-distribution parameters
    (N0, lambda, alpha) and evaluates :func:`~pyPIPS.DSDlib.calc_binned_DSD_from_params`
    at the Parsivel midpoint diameters.  The result is time-weighted and
    averaged into each ``sample_time`` interval.

    Parameters
    ----------
    transect_ds : xr.Dataset
        Output of :func:`find_transect_grid_intersections`.
        Must have coords ``all_time`` and ``sample_time``.
    ds : xr.Dataset
        Normalised model Dataset with canonical field names
        (``rhoa``, ``qr``, ``ntr``, ``zr``).
        Scalar fields are expected to have dimensions
        ``(time, [vertical,] yc, xc)``.
    grid_idx_ds : xr.Dataset
        Output of :func:`get_transect_grid_indices`.
        Variables: ``i_idx``, ``j_idx``, optionally ``time_idx``.
        Coord: ``all_time``.
    Dmax : float or None
        Maximum diameter (m).  Uses the full bin range when ``None``.
    level : int
        Vertical level index (0 = lowest model level).

    Returns
    -------
    xr.Dataset
        Variables:

        ``ND_model`` – (sample_time, diameter_bin) number distribution in
        m\ :sup:`-3` mm\ :sup:`-1`, time-weighted mean over each interval.

        Coord: ``sample_time``.
    """
    Dmax, Dmax_index = get_Dmax_index(Dr, Dmax)
    nbins = Dmax_index + 1

    # Parsivel midpoint diameters (mm) for the bins up to Dmax.
    D_mm = D[:nbins] * 1000.0  # m → mm

    sampling_times = transect_ds.coords["sample_time"].values
    all_times = transect_ds.coords["all_time"].values

    dt = all_times[1:] - all_times[:-1]          # durations between all_time points (s)
    sampling_dt = sampling_times[1:] - sampling_times[:-1]  # sample-interval durations (s)

    i_idx_arr = grid_idx_ds["i_idx"].values.astype(int)
    j_idx_arr = grid_idx_ds["j_idx"].values.astype(int)
    has_time_idx = "time_idx" in grid_idx_ds
    if has_time_idx:
        time_idx_arr = grid_idx_ds["time_idx"].values.astype(int)

    ntimes = len(all_times)
    time_coord = "time"

    # ------------------------------------------------------------------
    # Bulk numpy extraction of model scalars along the transect.
    # ------------------------------------------------------------------

    def _get_field_numpy(name: str) -> np.ndarray:
        vals = ds[name].values
        if vals.ndim == 4:          # (time, level, y, x)
            return vals[:, level, :, :]
        if vals.ndim == 3 and time_coord in ds[name].dims:
            return vals             # (time, y, x)
        if vals.ndim == 3:          # (level, y, x) — no time dimension
            return vals[level]      # → (y, x)
        return vals                 # (y, x)

    def _index_field(arr: np.ndarray, tidx) -> np.ndarray:
        if arr.ndim == 3:
            return arr[tidx, j_idx_arr, i_idx_arr]
        return arr[j_idx_arr, i_idx_arr]

    if has_time_idx:
        tidx = time_idx_arr
    elif time_coord in ds.dims:
        tidx = np.zeros(ntimes, dtype=int)
    else:
        tidx = None

    # TODO: leverage model_config variable mapping for field names.
    # EDIT: changed this to use the canonical field names in the model dataset, so it should work
    # assuming that the model dataset has been "normalized" to have the canonical field names
    # (rhoa, qr, ntr, zr).
    if "rhoa" in ds:
        rhoa_arr = _index_field(_get_field_numpy("rhoa"), tidx)
    elif "rho" in ds:
        rhoa_arr = _index_field(_get_field_numpy("rho"), tidx)
    else:
        p_arr  = _index_field(_get_field_numpy("p"), tidx)
        pt_arr = _index_field(_get_field_numpy("th"),  tidx)
        qv_arr = _index_field(_get_field_numpy("qv"),  tidx)
        rhoa_arr = np.array([
            float(thermo.calrho(p_arr[n], pt_arr[n], qv_arr[n]))
            for n in range(ntimes)
        ])

    qr_arr  = _index_field(_get_field_numpy("qr"),  tidx)
    ntr_arr = _index_field(_get_field_numpy("ntr"), tidx)
    zr_arr  = _index_field(_get_field_numpy("zr"), tidx)

    # ------------------------------------------------------------------
    # Gamma-distribution parameters from model moments.
    # ------------------------------------------------------------------
    alphar_arr = np.array([
        dualpol.solve_alpha_iter(
            rhoa_arr[n], mur, qr_arr[n], ntr_arr[n], zr_arr[n], rhorcst
        ).squeeze()
        for n in range(ntimes)
    ])

    N0r   = dsd.calc_N0_gamma(rhoa_arr, qr_arr, ntr_arr, cr, alphar_arr)
    lamdar = dsd.calc_lamda_gamma(rhoa_arr, qr_arr, ntr_arr, cr, alphar_arr)

    # ------------------------------------------------------------------
    # Evaluate the gamma DSD at all all_time points simultaneously.
    # Broadcasting: params are (ntimes, 1) × D_mm is (1, nbins).
    # calc_binned_DSD_from_params divides D by 1000 internally (mm → m),
    # and returns N(D) in m^-4.  Multiply by 1e-3 → m^-3 mm^-1.
    # ------------------------------------------------------------------
    N0r_np    = np.asarray(N0r)[:, None]       # (ntimes, 1)
    lamdar_np = np.asarray(lamdar)[:, None]    # (ntimes, 1)
    alpha_np  = alphar_arr[:, None]            # (ntimes, 1)
    D_bcast   = D_mm[None, :]                  # (1, nbins)

    ND_all_times = 1.0e-3 * dsd.calc_binned_DSD_from_params(
        N0r_np, lamdar_np, alpha_np, D_bcast
    )  # (ntimes, nbins), m^-3 mm^-1

    # Replace NaN/inf values (dry cells where qr≈0) with 0.
    ND_all_times = np.where(np.isfinite(ND_all_times), ND_all_times, 0.0)

    # ------------------------------------------------------------------
    # Time-weighted average of N(D) into each sample_time interval.
    # ------------------------------------------------------------------
    nsamples = len(sampling_times)
    ND_model = np.zeros((nsamples, nbins))

    sample_indices = np.searchsorted(all_times, sampling_times, side="left")

    for s, sample_index in enumerate(sample_indices[:-1]):
        sample_index_end = sample_indices[s + 1]
        current_slice = slice(sample_index, sample_index_end)
        weights = dt[current_slice]          # (n_sub,) seconds
        ND_slice = ND_all_times[current_slice]  # (n_sub, nbins)
        if sampling_dt[s] > 0:
            ND_model[s + 1, :] = (
                np.sum(ND_slice * weights[:, None], axis=0) / sampling_dt[s]
            )

    ND_model = np.ma.masked_invalid(ND_model)

    return xr.Dataset(
        {
            "ND_model": xr.DataArray(
                ND_model, dims=["sample_time", "diameter_bin"]
            ),
        },
        coords={"sample_time": sampling_times},
    )


def calc_model_PSD_along_transect_hail(
    transect_ds: xr.Dataset,
    ds: xr.Dataset,
    grid_idx_ds: xr.Dataset,
    Dmax: float | None = None,
    level: int = 0,
) -> xr.Dataset:
    """
    Compute the raw model gamma PSD along a transect, binned into Parsivel bins.

    Unlike :func:`sample_model_PSD_along_transect`, this function does not run
    the Parsivel simulator.  Instead it uses the predicted model moments
    (``qhl``, ``nthl``, ``zhl``) to derive the gamma-distribution parameters
    (N0, lambda, alpha) and evaluates :func:`~pyPIPS.DSDlib.calc_binned_DSD_from_params`
    at the Parsivel midpoint diameters.  The result is time-weighted and
    averaged into each ``sample_time`` interval.

    Parameters
    ----------
    transect_ds : xr.Dataset
        Output of :func:`find_transect_grid_intersections`.
        Must have coords ``all_time`` and ``sample_time``.
    ds : xr.Dataset
        Normalised model Dataset with canonical field names
        (``rhoa``, ``qhl``, ``nthl``, ``zhl``).
        Scalar fields are expected to have dimensions
        ``(time, [vertical,] yc, xc)``.
    grid_idx_ds : xr.Dataset
        Output of :func:`get_transect_grid_indices`.
        Variables: ``i_idx``, ``j_idx``, optionally ``time_idx``.
        Coord: ``all_time``.
    Dmax : float or None
        Maximum diameter (m).  Uses the full bin range when ``None``.
    level : int
        Vertical level index (0 = lowest model level).

    Returns
    -------
    xr.Dataset
        Variables:

        ``ND_model`` – (sample_time, diameter_bin) number distribution in
        m\ :sup:`-3` mm\ :sup:`-1`, time-weighted mean over each interval.

        Coord: ``sample_time``.
    """
    Dmax, Dmax_index = get_Dmax_index(Dr, Dmax)
    nbins = Dmax_index + 1

    # Parsivel midpoint diameters (mm) for the bins up to Dmax.
    D_mm = D[:nbins] * 1000.0  # m → mm

    sampling_times = transect_ds.coords["sample_time"].values
    all_times = transect_ds.coords["all_time"].values

    dt = all_times[1:] - all_times[:-1]          # durations between all_time points (s)
    sampling_dt = sampling_times[1:] - sampling_times[:-1]  # sample-interval durations (s)

    i_idx_arr = grid_idx_ds["i_idx"].values.astype(int)
    j_idx_arr = grid_idx_ds["j_idx"].values.astype(int)
    has_time_idx = "time_idx" in grid_idx_ds
    if has_time_idx:
        time_idx_arr = grid_idx_ds["time_idx"].values.astype(int)

    ntimes = len(all_times)
    time_coord = "time"

    #------------------------------------------------------------------
    # Bulk numpy extraction of model scalars along the transect.
    #------------------------------------------------------------------

    def _get_field_numpy(name: str) -> np.ndarray:
        vals = ds[name].values
        if vals.ndim == 4:          # (time, level, y, x)
            return vals[:, level, :, :]
        if vals.ndim == 3 and time_coord in ds[name].dims:
            return vals             # (time, y, x)
        if vals.ndim == 3:          # (level, y, x) — no time dimension
            return vals[level]      # → (y, x)
        return vals                 # (y, x)

    def _index_field(arr: np.ndarray, tidx) -> np.ndarray:
        if arr.ndim == 3:
            return arr[tidx, j_idx_arr, i_idx_arr]
        return arr[j_idx_arr, i_idx_arr]

    if has_time_idx:
        tidx = time_idx_arr
    elif time_coord in ds.dims:
        tidx = np.zeros(ntimes, dtype=int)
    else:
        tidx = None

    # TODO: leverage model_config variable mapping for field names.
    if "rhoa" in ds:
        rhoa_arr = _index_field(_get_field_numpy("rhoa"), tidx)
    elif "rho" in ds:
        rhoa_arr = _index_field(_get_field_numpy("rho"), tidx)
    else:
        p_arr  = _index_field(_get_field_numpy("p"), tidx)
        pt_arr = _index_field(_get_field_numpy("th"),  tidx)
        qv_arr = _index_field(_get_field_numpy("qv"),  tidx)
        rhoa_arr = np.array([
            float(thermo.calrho(p_arr[n], pt_arr[n], qv_arr[n]))
            for n in range(ntimes)
        ])

    qhl_arr  = _index_field(_get_field_numpy("qh"), tidx)
    nthl_arr = _index_field(_get_field_numpy("nth"), tidx)
    zhl_arr  = _index_field(_get_field_numpy("zh"), tidx)
    vhl_arr  = _index_field(_get_field_numpy("vh"), tidx)

    # calculate density of hail
    rhohl_arr = qhl_arr/vhl_arr
    chl_arr = np.pi / 6. * rhohl_arr

    # ------------------------------------------------------------------
    # Derive alphahl from reflectivity using the iterative solver.
    # ------------------------------------------------------------------
    alphahl_arr = np.array([
        dualpol.solve_alpha_iter(
            rhoa_arr[n], muhl, qhl_arr[n], nthl_arr[n], zhl_arr[n], rhohl_arr[n]
        ).squeeze()
        for n in range(ntimes)
    ])
    # ------------------------------------------------------------------
    # Derived DSD parameters.
    # ------------------------------------------------------------------
    N0hl   = dsd.calc_N0_gamma(rhoa_arr, qhl_arr, nthl_arr, chl_arr, alphahl_arr)
    lamdahl = dsd.calc_lamda_gamma(rhoa_arr, qhl_arr, nthl_arr, chl_arr, alphahl_arr)
    # ------------------------------------------------------------------
    # Evaluate the gamma DSD at all all_time points simultaneously.
    # Broadcasting: params are (ntimes, 1) × D_mm is (1, nbins).
    # calc_binned_DSD_from_params divides D by 1000 internally (mm → m),
    # and returns N(D) in m^-4.  Multiply by 1e-3 → m^-3 mm^-1.
    # ------------------------------------------------------------------
    N0hl_np    = np.asarray(N0hl)[:, None]       # (ntimes, 1)
    lamdahl_np = np.asarray(lamdahl)[:, None]    # (ntimes, 1)
    alpha_np  = alphahl_arr[:, None]            # (ntimes, 1)
    D_bcast   = D_mm[None, :]                  # (1, nbins)

    ND_all_times = 1.0e-3 * dsd.calc_binned_DSD_from_params(
        N0hl_np, lamdahl_np, alpha_np, D_bcast
    )  # (ntimes, nbins), m^-3 mm^-1

    # Replace NaN/inf values (dry cells where qhl≈0) with 0.
    ND_all_times = np.where(np.isfinite(ND_all_times), ND_all_times, 0.0)
    print("ND_all_times:", ND_all_times)

    # ------------------------------------------------------------------
    # Time-weighted average of N(D) into each sample_time interval.
    # ------------------------------------------------------------------
    nsamples = len(sampling_times)
    ND_model = np.zeros((nsamples, nbins))

    sample_indices = np.searchsorted(all_times, sampling_times, side="left")

    for s, sample_index in enumerate(sample_indices[:-1]):
        sample_index_end = sample_indices[s + 1]
        current_slice = slice(sample_index, sample_index_end)
        weights = dt[current_slice]          # (n_sub,) seconds
        ND_slice = ND_all_times[current_slice]  # (n_sub, nbins)
        if sampling_dt[s] > 0:
            ND_model[s + 1, :] = (
                np.sum(ND_slice * weights[:, None], axis=0) / sampling_dt[s]
            )

    ND_model = np.ma.masked_invalid(ND_model)

    return xr.Dataset(
        {
            "ND_model": xr.DataArray(
                ND_model, dims=["sample_time", "diameter_bin"]
            ),
        },
        coords={"sample_time": sampling_times},
    )



def sample_model_PSD_along_transect_pressure(
    transect_ds: xr.Dataset,
    ds: xr.Dataset,
    grid_idx_ds: xr.Dataset,
    Dmax: float | None = None,
    level: int = 0,
) -> xr.Dataset:
    """
    Sample the model pressure along a transect using the Parsivel simulator.
    """

    Dmax, Dmax_index = get_Dmax_index(Dr, Dmax)

    sampling_times = transect_ds.coords["sample_time"].values
    all_times = transect_ds.coords["all_time"].values

    # Time deltas between consecutive all_time points (seconds).
    dt = all_times[1:] - all_times[:-1]
    # Per-sample-interval durations (seconds), implicit in sampling_times.
    sampling_dt = sampling_times[1:] - sampling_times[:-1]

    i_idx_arr = grid_idx_ds["i_idx"].values.astype(int)
    j_idx_arr = grid_idx_ds["j_idx"].values.astype(int)
    has_time_idx = "time_idx" in grid_idx_ds
    if has_time_idx:
        time_idx_arr = grid_idx_ds["time_idx"].values.astype(int)

    ntimes = len(all_times)
    time_coord = "time"

    # ------------------------------------------------------------------
    # Extract model scalars for every all_time point using bulk numpy
    # indexing to avoid per-step xarray overhead.
    # ------------------------------------------------------------------

    def _get_field_numpy(name: str) -> np.ndarray:
        """Return a 3-D (time, y, x) or 2-D (y, x) numpy array for *name*."""
        vals = ds[name].values
        if vals.ndim == 4:          # (time, level, y, x)
            return vals[:, level, :, :]
        if vals.ndim == 3 and time_coord in ds[name].dims:
            return vals             # (time, y, x)
        if vals.ndim == 3:          # (level, y, x) — no time dimension
            return vals[level]      # → (y, x)
        return vals                 # (y, x)

    def _index_field(arr: np.ndarray, tidx: np.ndarray | int) -> np.ndarray:
        """Advanced-index *arr* to extract the ntimes-long 1-D transect."""
        if arr.ndim == 3:
            return arr[tidx, j_idx_arr, i_idx_arr]
        # 2-D (y, x) — same spatial location for all times
        return arr[j_idx_arr, i_idx_arr]

    if has_time_idx:
        tidx = time_idx_arr
    elif time_coord in ds.dims:
        tidx = np.zeros(ntimes, dtype=int)
    else:
        tidx = None  # no time dimension

    p_arr  = _index_field(_get_field_numpy("prs"), tidx)


    return xr.Dataset(
        {
            "p": xr.DataArray(
                p_arr, dims=["all_time"]
            ),
        },
        coords={"all_time": all_times},
    )

def sample_model_PSD_along_transect_theta(
    transect_ds: xr.Dataset,
    ds: xr.Dataset,
    grid_idx_ds: xr.Dataset,
    Dmax: float | None = None,
    level: int = 0,
) -> xr.Dataset:
    """
    Sample the model potential temperature along a transect using the Parsivel simulator.
    """

    Dmax, Dmax_index = get_Dmax_index(Dr, Dmax)

    sampling_times = transect_ds.coords["sample_time"].values
    all_times = transect_ds.coords["all_time"].values

    # Time deltas between consecutive all_time points (seconds).
    dt = all_times[1:] - all_times[:-1]
    # Per-sample-interval durations (seconds), implicit in sampling_times.
    sampling_dt = sampling_times[1:] - sampling_times[:-1]

    i_idx_arr = grid_idx_ds["i_idx"].values.astype(int)
    j_idx_arr = grid_idx_ds["j_idx"].values.astype(int)
    has_time_idx = "time_idx" in grid_idx_ds
    if has_time_idx:
        time_idx_arr = grid_idx_ds["time_idx"].values.astype(int)

    ntimes = len(all_times)
    time_coord = "time"

    # ------------------------------------------------------------------
    # Extract model scalars for every all_time point using bulk numpy
    # indexing to avoid per-step xarray overhead.
    # ------------------------------------------------------------------

    def _get_field_numpy(name: str) -> np.ndarray:
        """Return a 3-D (time, y, x) or 2-D (y, x) numpy array for *name*."""
        vals = ds[name].values
        if vals.ndim == 4:          # (time, level, y, x)
            return vals[:, level, :, :]
        if vals.ndim == 3 and time_coord in ds[name].dims:
            return vals             # (time, y, x)
        if vals.ndim == 3:          # (level, y, x) — no time dimension
            return vals[level]      # → (y, x)
        return vals                 # (y, x)

    def _index_field(arr: np.ndarray, tidx: np.ndarray | int) -> np.ndarray:
        """Advanced-index *arr* to extract the ntimes-long 1-D transect."""
        if arr.ndim == 3:
            return arr[tidx, j_idx_arr, i_idx_arr]
        # 2-D (y, x) — same spatial location for all times
        return arr[j_idx_arr, i_idx_arr]

    if has_time_idx:
        tidx = time_idx_arr
    elif time_coord in ds.dims:
        tidx = np.zeros(ntimes, dtype=int)
    else:
        tidx = None  # no time dimension

    th_arr  = _index_field(_get_field_numpy("th"), tidx)


    return xr.Dataset(
        {
            "th": xr.DataArray(
                th_arr, dims=["all_time"]
            ),
        },
        coords={"all_time": all_times},
    )

def sample_model_PSD_along_transect_wind(
    transect_ds: xr.Dataset,
    ds: xr.Dataset,
    grid_idx_ds: xr.Dataset,
    Dmax: float | None = None,
    level: int = 0,
) -> xr.Dataset:
    """
    Sample the model surface winds along a transect using the Parsivel simulator.
    """

    Dmax, Dmax_index = get_Dmax_index(Dr, Dmax)

    sampling_times = transect_ds.coords["sample_time"].values
    all_times = transect_ds.coords["all_time"].values

    # Time deltas between consecutive all_time points (seconds).
    dt = all_times[1:] - all_times[:-1]
    # Per-sample-interval durations (seconds), implicit in sampling_times.
    sampling_dt = sampling_times[1:] - sampling_times[:-1]

    i_idx_arr = grid_idx_ds["i_idx"].values.astype(int)
    j_idx_arr = grid_idx_ds["j_idx"].values.astype(int)
    has_time_idx = "time_idx" in grid_idx_ds
    if has_time_idx:
        time_idx_arr = grid_idx_ds["time_idx"].values.astype(int)

    ntimes = len(all_times)
    time_coord = "time"

    # ------------------------------------------------------------------
    # Extract model scalars for every all_time point using bulk numpy
    # indexing to avoid per-step xarray overhead.
    # ------------------------------------------------------------------

    def _get_field_numpy(name: str) -> np.ndarray:
        """Return a 3-D (time, y, x) or 2-D (y, x) numpy array for *name*."""
        vals = ds[name].values
        if vals.ndim == 4:          # (time, level, y, x)
            return vals[:, level, :, :]
        if vals.ndim == 3 and time_coord in ds[name].dims:
            return vals             # (time, y, x)
        if vals.ndim == 3:          # (level, y, x) — no time dimension
            return vals[level]      # → (y, x)
        return vals                 # (y, x)

    def _index_field(arr: np.ndarray, tidx: np.ndarray | int) -> np.ndarray:
        """Advanced-index *arr* to extract the ntimes-long 1-D transect."""
        if arr.ndim == 3:
            return arr[tidx, j_idx_arr, i_idx_arr]
        # 2-D (y, x) — same spatial location for all times
        return arr[j_idx_arr, i_idx_arr]

    if has_time_idx:
        tidx = time_idx_arr
    elif time_coord in ds.dims:
        tidx = np.zeros(ntimes, dtype=int)
    else:
        tidx = None  # no time dimension

    u  = _index_field(_get_field_numpy("u"), tidx)
    v  = _index_field(_get_field_numpy("v"), tidx)


    return xr.Dataset(
        {
            "u": xr.DataArray(
                u, dims=["all_time"]
            ),
            "v": xr.DataArray(
                v, dims=["all_time"]
            ),
        },
        coords={"all_time": all_times},
    )

def sample_model_PSD_along_transect_qv(
    transect_ds: xr.Dataset,
    ds: xr.Dataset,
    grid_idx_ds: xr.Dataset,
    Dmax: float | None = None,
    level: int = 0,
) -> xr.Dataset:
    """
    Sample the model vapor mixing ratio along a transect using the Parsivel simulator.
    """

    Dmax, Dmax_index = get_Dmax_index(Dr, Dmax)

    sampling_times = transect_ds.coords["sample_time"].values
    all_times = transect_ds.coords["all_time"].values

    # Time deltas between consecutive all_time points (seconds).
    dt = all_times[1:] - all_times[:-1]
    # Per-sample-interval durations (seconds), implicit in sampling_times.
    sampling_dt = sampling_times[1:] - sampling_times[:-1]

    i_idx_arr = grid_idx_ds["i_idx"].values.astype(int)
    j_idx_arr = grid_idx_ds["j_idx"].values.astype(int)
    has_time_idx = "time_idx" in grid_idx_ds
    if has_time_idx:
        time_idx_arr = grid_idx_ds["time_idx"].values.astype(int)

    ntimes = len(all_times)
    time_coord = "time"

    # ------------------------------------------------------------------
    # Extract model scalars for every all_time point using bulk numpy
    # indexing to avoid per-step xarray overhead.
    # ------------------------------------------------------------------

    def _get_field_numpy(name: str) -> np.ndarray:
        """Return a 3-D (time, y, x) or 2-D (y, x) numpy array for *name*."""
        vals = ds[name].values
        if vals.ndim == 4:          # (time, level, y, x)
            return vals[:, level, :, :]
        if vals.ndim == 3 and time_coord in ds[name].dims:
            return vals             # (time, y, x)
        if vals.ndim == 3:          # (level, y, x) — no time dimension
            return vals[level]      # → (y, x)
        return vals                 # (y, x)

    def _index_field(arr: np.ndarray, tidx: np.ndarray | int) -> np.ndarray:
        """Advanced-index *arr* to extract the ntimes-long 1-D transect."""
        if arr.ndim == 3:
            return arr[tidx, j_idx_arr, i_idx_arr]
        # 2-D (y, x) — same spatial location for all times
        return arr[j_idx_arr, i_idx_arr]

    if has_time_idx:
        tidx = time_idx_arr
    elif time_coord in ds.dims:
        tidx = np.zeros(ntimes, dtype=int)
    else:
        tidx = None  # no time dimension

    qv  = _index_field(_get_field_numpy("qv"), tidx)


    return xr.Dataset(
        {
            "qv": xr.DataArray(
                qv, dims=["all_time"]
            ),
        },
        coords={"all_time": all_times},
    )