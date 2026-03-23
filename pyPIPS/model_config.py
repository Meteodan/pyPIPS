"""Model-specific configuration and variable-name mappings.

All downstream simulator functions work with a *canonical* set of variable and
coordinate names.  A :class:`ModelConfig` provides the mapping from model-native
names to these canonical names, and :func:`open_and_normalize` applies that
mapping to a raw ``xr.Dataset`` at load time so that every downstream function
can use, e.g., ``ds['qr']`` or ``ds.coords['xc']`` regardless of which model
produced the file.

Canonical coordinate names
---------------------------
xc, yc, zc   : cell-center coordinates (m)
xe, ye, ze   : cell-edge coordinates (m)
time         : time coordinate (s from simulation start, or numeric time unit
               stored in the file - not decoded to datetime64 by default)

Canonical field names (representative subset)
----------------------------------------------
U, V, W      : wind components (staggered or destaggered, model-dependent)
TH           : potential temperature (K)
QV           : water-vapour mixing ratio (kg kg-1)
P            : total pressure (Pa)
qr           : rain mixing ratio (kg kg-1)
ntr          : rain total number concentration (m-3)
zr           : rain reflectivity moment (m6 m-3)
qh, nth      : hail mixing ratio / number concentration
qg, ntg      : graupel mixing ratio / number concentration
rhoa         : air density (kg m-3)
alphar, alphah, alphag : shape parameters for rain/hail/graupel
DBZ          : simulated equivalent reflectivity (dBZ)
"""

from __future__ import annotations

import glob as _glob
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
import xarray as xr


# ---------------------------------------------------------------------------
# Configuration dataclass
# ---------------------------------------------------------------------------

@dataclass
class ModelConfig:
    """
    Holds all information needed to load and normalize a model dataset.

    Parameters
    ----------
    model_type : str
        Human-readable label, e.g. ``'CM1'``, ``'COMMAS'``, ``'WRF'``.
    coord_map : dict
        Canonical-name → native-coordinate-name mapping.
    var_map : dict
        Canonical-name → native-variable-name mapping.
    file_pattern : str or list of str
        Glob pattern(s) or explicit file paths for :func:`xr.open_mfdataset`.
    engine : str
        NetCDF engine (default ``'netcdf4'``).
    chunks : dict or None
        Dask chunking specification, e.g. ``{'time': 1}``.
    combine : str
        How to combine multiple files (``'by_coords'`` or ``'nested'``).
    decode_times : bool
        Whether xarray should decode times.  Set ``False`` (default) to keep
        numeric seconds as stored by CM1/COMMAS.
    preprocess : callable or None
        Optional pre-processing function applied to each file before merging;
        passed directly to :func:`xr.open_mfdataset`.
    run_name : str
        Optional run label stored in ``ds.attrs['run_name']``.
    coord_to_meters_factors : dict or None
        Explicit conversion factors for spatial coordinates to meters, keyed by
        canonical name (e.g. ``{'xc': 0.001, 'yc': 0.001}`` to convert from km).
        If not provided, units attributes on coordinates are checked.
    """

    model_type: str
    coord_map: dict[str, str] = field(default_factory=dict)
    var_map: dict[str, str] = field(default_factory=dict)
    file_pattern: str | list[str] = ""
    engine: str = "netcdf4"
    chunks: dict[str, int] | None = None
    combine: str = "by_coords"
    decode_times: bool = False
    preprocess: Callable[[xr.Dataset], xr.Dataset] | None = None
    run_name: str = ""
    coord_to_meters_factors: dict[str, float] | None = None


# ---------------------------------------------------------------------------
# Pre-built coordinate / variable maps for supported models
# ---------------------------------------------------------------------------

#: CM1 netCDF output (default namelist.input settings, cm1out_*.nc)
CM1_COORD_MAP: dict[str, str] = {
    "xc": "xh",   # scalar (cell-center) x in meters
    "yc": "yh",
    "zc": "zh",
    "xe": "xf",   # cell-edge (staggered) x
    "ye": "yf",
    "ze": "zf",
    "time": "time",
}

CM1_VAR_MAP: dict[str, str] = {
    "u": "uinterp",   # wind interpolated to scalar points when available
    "v": "vinterp",
    "w": "winterp",
    "th": "th",
    "qv": "qv",
    "p": "prs",       # total pressure in Pa
    "qr": "qr",
    "ntr": "nrain",
    "zr": "zrain",
    "qh": "qhl",      # hail-liquid; name varies with microphysics scheme
    "nth": "nhl",
    "qg": "qice",     # approximate; adjust for your scheme
    "ntg": "nice",
    "rhoa": "rho",
    "DBZ": "dbz",
}

#: COMMAS netCDF output
COMMAS_COORD_MAP: dict[str, str] = {
    "xc": "xc",
    "yc": "yc",
    "zc": "zc",
    "xe": "xe",
    "ye": "ye",
    "ze": "ze",
    "time": "time",
}

COMMAS_VAR_MAP: dict[str, str] = {
    "u": "U",
    "v": "V",
    "w": "W",
    "th": "TH",
    "qv": "QV",
    "p": "P",
    "qr": "qr",
    "ntr": "ntr",
    "zr": "zr",
    "qh": "qh",
    "nth": "nth",
    "qg": "qg",
    "ntg": "ntg",
    "rhoa": "rhoa",
    "alphar": "alphar",
    "alphah": "alphah",
    "alphag": "alphag",
    "DBZ": "DBZ",
}

#: WRF (placeholder - extend as needed for your configuration)
WRF_COORD_MAP: dict[str, str] = {
    "xc": "west_east",
    "yc": "south_north",
    "zc": "bottom_top",
    "time": "Time",
}

WRF_VAR_MAP: dict[str, str] = {
    "u": "U",
    "v": "V",
    "w": "W",
    "th": "T",          # perturbation θ; add 300 K downstream
    "qv": "QVAPOR",
    "p": "P",           # perturbation pressure; add PB downstream
    "qr": "QRAIN",
    "DBZ": "REFL_10CM",
}

_MODEL_REGISTRY: dict[str, tuple[dict, dict]] = {
    "CM1": (CM1_COORD_MAP, CM1_VAR_MAP),
    "COMMAS": (COMMAS_COORD_MAP, COMMAS_VAR_MAP),
    "WRF": (WRF_COORD_MAP, WRF_VAR_MAP),
}


def make_config(
    model_type: str,
    file_pattern: str | list[str] = "",
    *,
    extra_coord_map: dict[str, str] | None = None,
    extra_var_map: dict[str, str] | None = None,
    engine: str = "netcdf4",
    chunks: dict[str, int] | None = None,
    combine: str = "by_coords",
    decode_times: bool = False,
    run_name: str = "",
) -> ModelConfig:
    """
    Convenience constructor that looks up built-in coord/var maps for
    *model_type* and merges any caller-supplied overrides.

    Parameters
    ----------
    model_type : {'CM1', 'COMMAS', 'WRF'}
    file_pattern : str or list of str, optional
        Glob pattern(s) used by :func:`open_and_normalize`.  Not required
        when you already have a dataset and only need :func:`normalize`.
    extra_coord_map : dict, optional
        Override or extend the built-in coordinate map.
    extra_var_map : dict, optional
        Override or extend the built-in variable map.
    engine, chunks, combine, decode_times, run_name
        Forwarded to :class:`ModelConfig`.

    Returns
    -------
    ModelConfig
    """
    if model_type not in _MODEL_REGISTRY:
        msg = (
            f"Unknown model_type {model_type!r}. "
            f"Choose from {list(_MODEL_REGISTRY)} or build a ModelConfig manually."
        )
        raise ValueError(msg)
    coord_map, var_map = _MODEL_REGISTRY[model_type]
    coord_map = dict(coord_map)
    var_map = dict(var_map)
    if extra_coord_map:
        coord_map.update(extra_coord_map)
    if extra_var_map:
        var_map.update(extra_var_map)
    return ModelConfig(
        model_type=model_type,
        coord_map=coord_map,
        var_map=var_map,
        file_pattern=file_pattern,
        engine=engine,
        chunks=chunks,
        combine=combine,
        decode_times=decode_times,
        run_name=run_name,
    )


# ---------------------------------------------------------------------------
# Dataset loading and normalisation
# ---------------------------------------------------------------------------

def _invert_map(m: dict[str, str]) -> dict[str, str]:
    """Return ``{native: canonical}`` from ``{canonical: native}``."""
    return {v: k for k, v in m.items()}


def _get_unit_to_meters_factor(unit_str: str) -> float | None:
    """
    Infer a conversion factor to meters from a units string.

    Parameters
    ----------
    unit_str : str
        Units string, e.g. 'km', 'meters', 'm', 'centimeters'.

    Returns
    -------
    float or None
        Conversion factor to multiply by to get meters, or None if not recognized.
    """
    unit_lower = unit_str.lower().strip()
    # Normalize common abbreviations and variations
    if unit_lower in ('km', 'kilometer', 'kilometers'):
        return 1000.0
    if unit_lower in ('m', 'meter', 'meters'):
        return 1.0
    if unit_lower in ('cm', 'centimeter', 'centimeters'):
        return 0.01
    if unit_lower in ('mm', 'millimeter', 'millimeters'):
        return 0.001
    return None


def _convert_spatial_coords_to_meters(ds: xr.Dataset, cfg: ModelConfig) -> xr.Dataset:
    """
    Convert spatial coordinate arrays to meters.

    First checks explicit conversion factors in ``cfg.coord_to_meters_factors``.
    If not found, reads the ``units`` attribute of each coordinate and infers
    the conversion factor.

    Preserves all existing coordinate attributes, updating only the 'units'
    attribute to 'meters' when a conversion is applied.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with already-renamed (canonical) coordinates.
    cfg : ModelConfig

    Returns
    -------
    xr.Dataset
        Dataset with coordinates converted to meters (attributes preserved).
    """
    explicit_factors = cfg.coord_to_meters_factors or {}
    spatial_coords = ('xc', 'yc', 'zc', 'xe', 'ye', 'ze')
    coords_to_update = {}

    for coord_name in spatial_coords:
        if coord_name not in ds.coords:
            continue

        # Determine conversion factor: explicit > inferred from units > no conversion
        factor = explicit_factors.get(coord_name)
        if factor is None:
            # Try to infer from units attribute
            units_attr = ds.coords[coord_name].attrs.get('units')
            if units_attr:
                factor = _get_unit_to_meters_factor(units_attr)

        # Apply conversion if factor != 1.0
        if factor is not None and factor != 1.0:
            # Get original coordinate and preserve all attributes
            orig_coord = ds.coords[coord_name]
            coord_values = orig_coord.values * factor

            # Create new coordinate with updated values and preserved attributes
            new_attrs = dict(orig_coord.attrs)
            new_attrs['units'] = 'meters'
            new_coord = xr.DataArray(coord_values, dims=orig_coord.dims, attrs=new_attrs)
            coords_to_update[coord_name] = new_coord

    # Apply all coordinate updates at once using assign_coords
    if coords_to_update:
        ds = ds.assign_coords(coords_to_update)

    return ds


def normalize(
    ds: xr.Dataset,
    cfg: ModelConfig,
    *,
    rename_vars: bool = True,
    rename_coords: bool = True,
) -> xr.Dataset:
    """
    normalize an already-loaded model dataset to canonical names.

    Use this when you have already opened the dataset yourself (e.g. with
    ``xr.open_dataset`` or ``xr.open_mfdataset``) and only need the
    rename/metadata step.

    After this call every downstream function can use, e.g., ``ds['qr']``,
    ``ds.coords['xc']`` regardless of which model produced the file.

    Grid-spacing attributes ``dx``, ``dy`` and ``dz`` (where uniform) are
    computed from the coordinate arrays and stored in ``ds.attrs``.

    Parameters
    ----------
    ds : xr.Dataset
        Raw model dataset with native variable/coordinate names.
    cfg : ModelConfig
        Configuration supplying the coordinate and variable name maps.
        The ``file_pattern``, ``engine``, ``chunks``, ``combine``,
        ``decode_times``, and ``preprocess`` fields are ignored.
    rename_vars : bool, optional
        Apply ``cfg.var_map`` to rename data variables.  Default ``True``.
    rename_coords : bool, optional
        Apply ``cfg.coord_map`` to rename coordinates.  Default ``True``.

    Returns
    -------
    xr.Dataset
        Dataset with canonical names and grid-spacing attrs.
    """
    # Build rename map: native → canonical, but only for names present in ds.
    combined: dict[str, str] = {}
    if rename_coords:
        combined.update(cfg.coord_map)
    if rename_vars:
        combined.update(cfg.var_map)
    full_native_to_canonical = _invert_map(combined)
    rename_map = {
        native: canonical
        for native, canonical in full_native_to_canonical.items()
        if native in ds or native in ds.coords
    }
    if rename_map:
        ds = ds.rename(rename_map)

    # Convert spatial coordinates to meters.
    ds = _convert_spatial_coords_to_meters(ds, cfg)

    # Store configuration metadata as attrs.
    ds.attrs["model_type"] = cfg.model_type
    if cfg.run_name:
        ds.attrs["run_name"] = cfg.run_name

    # Compute and cache grid spacings.
    for axis, coord_name in (("dx", "xc"), ("dy", "yc"), ("dz", "zc")):
        if coord_name in ds.coords:
            vals = np.asarray(ds.coords[coord_name].values).ravel()
            if vals.size > 1:
                ds.attrs[axis] = float(np.median(np.diff(vals)))

    return ds


def open_and_normalize(
    cfg: ModelConfig,
    *,
    rename_vars: bool = True,
    rename_coords: bool = True,
) -> xr.Dataset:
    """
    Open the model files described by *cfg* and rename all coordinates and
    variables to the canonical names defined in ``cfg.coord_map`` /
    ``cfg.var_map``.

    After this call every downstream function can use, e.g., ``ds['qr']``,
    ``ds.coords['xc']`` regardless of which model produced the file.

    Grid-spacing attributes ``dx``, ``dy`` and ``dz`` (where uniform) are
    computed from the coordinate arrays and stored in ``ds.attrs``.

    Parameters
    ----------
    cfg : ModelConfig
    rename_vars : bool, optional
        Apply ``cfg.var_map`` to rename data variables.  Default ``True``.
    rename_coords : bool, optional
        Apply ``cfg.coord_map`` to rename coordinates.  Default ``True``.

    Returns
    -------
    xr.Dataset
        Dataset with canonical names and grid-spacing attrs.

    Raises
    ------
    FileNotFoundError
        If no files are found matching ``cfg.file_pattern``.
    """
    pattern = cfg.file_pattern
    if isinstance(pattern, str):
        paths = sorted(_glob.glob(pattern))
        if not paths:
            msg = f"No model files matched pattern: {pattern!r}"
            raise FileNotFoundError(msg)
    else:
        paths = list(pattern)

    ds = xr.open_mfdataset(
        paths,
        engine=cfg.engine,
        chunks=cfg.chunks,
        combine=cfg.combine,
        decode_times=cfg.decode_times,
        preprocess=cfg.preprocess,
    )

    return normalize(ds, cfg, rename_vars=rename_vars, rename_coords=rename_coords)


# ---------------------------------------------------------------------------
# Grid array extraction
# ---------------------------------------------------------------------------

def get_grid_arrays(ds: xr.Dataset) -> dict[str, np.ndarray]:
    """
    Extract 1-D coordinate arrays from a normalized dataset.

    Cell-edge arrays are synthesised when not present (by assuming uniform
    spacing and placing edges at the midpoints between cell centers).

    Parameters
    ----------
    ds : xr.Dataset
        normalized model dataset (output of :func:`open_and_normalize`).

    Returns
    -------
    dict with keys:
        ``xc1d``, ``yc1d``, ``zc1d``   : cell-center coordinate vectors
        ``xe1d``, ``ye1d``, ``ze1d``   : cell-edge coordinate vectors
        ``dx``, ``dy``                 : median grid spacings (m)
    """

    def _to1d(name: str) -> np.ndarray | None:
        if name not in ds.coords:
            return None
        return np.asarray(ds.coords[name].values).ravel()

    def _edge(centers: np.ndarray) -> np.ndarray:
        """Build a len(centers)+1 edge array from cell centers."""
        d = np.diff(centers)
        left = centers[0] - 0.5 * d[0]
        right = centers[-1] + 0.5 * d[-1]
        return np.concatenate(([left], 0.5 * (centers[:-1] + centers[1:]), [right]))

    xc1d = _to1d("xc")
    yc1d = _to1d("yc")
    zc1d = _to1d("zc") if "zc" in ds.coords else np.array([0.0])

    xe1d = _to1d("xe") if "xe" in ds.coords else _edge(xc1d)
    ye1d = _to1d("ye") if "ye" in ds.coords else _edge(yc1d)
    ze1d = _to1d("ze") if "ze" in ds.coords else _edge(zc1d)

    dx = float(ds.attrs.get("dx", np.median(np.diff(xc1d))))
    dy = float(ds.attrs.get("dy", np.median(np.diff(yc1d))))

    return {
        "xc1d": xc1d,
        "yc1d": yc1d,
        "zc1d": zc1d,
        "xe1d": xe1d,
        "ye1d": ye1d,
        "ze1d": ze1d,
        "dx": dx,
        "dy": dy,
    }
