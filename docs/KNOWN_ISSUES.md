# Known Issues (for future review)

Issues found during code inventory that have been documented but intentionally not fixed yet.
Line numbers refer to the code as of October 4, 2026 and may drift.

## PIPS processing workflow

### ND and its QC variants use the original (pre-thermodynamic-QC) air density
- `ND`, `ND_qc`, `ND_roqc`, etc. are computed with a fall-speed density correction
  (`pips.calc_fallspeed_spectrum(..., correct_rho=True, rho=parsivel_combined_ds['rho'])`).
  The `rho` used is the value from `PIPS_csv_to_nc.py`, before any thermodynamic QC.
- The QC applied to the DSD (tags such as `qc`, `roqc`) is orthogonal to the QC applied to the
  conventional thermodynamic and kinematic data. As of October 2026, `apply_QC.py` recomputes
  `pt`, `qv`, and `rho` after thermodynamic QC, but `ND` and its variants are not recomputed with
  the updated `rho`. Fixing this properly needs a workflow-wide change (thermodynamic QC
  before ND construction, or recomputing ND afterward).
- Expected impact is small in almost all cases.

### `radar_to_PIPS.py`
- `--input-tag` (around lines 97-101): when set, filenames are never appended to
  `parsivel_combined_filenames`, so the list stays empty and nothing is processed.
- `--ngates2avg` is parsed but never passed to `interp_sweeps_to_one_PIPS_new`.
- The original radar sweep times are lost after `interp_like` onto the PIPS times, and the
  `elevation_angle` attribute records the requested `el_req`, not the actual elevation.

### `radarmodule.py`
- `interp_sweeps_to_one_PIPS_new` ignores its `compute_kdp` argument (hardcoded
  `compute_kdp=False` in the call to `readCFRadial_pyART`).
- `interp_one_sweep_to_one_PIPS`: when the nearest gate has `theta_index=0`, the
  neighborhood slice `theta_index-1:theta_index+2` is empty, so the gate average comes back empty.
- Beam height stored at the PIPS is the nearest gate center above the PIPS; beam width and
  vertical extent are not recorded.
- Some uncommitted changes to this module were lost earlier (see commit `e16fe46`).

### `calc_radar_retrieval_PIPS.py`
- Possible `KeyError`: the retrieval is called with `retrieval_tag=args.retrieval_tag`, but
  output keys are later looked up with a `retrieval_tag` that includes `QC_tag` unless the tag
  is `C08`/`Z01` or `--coefficients` was given.
- The radar-retrieved binned `ND_retr` is computed but discarded.

### Retrieval lookup tables
- `create_retrieval_lookup.py` writes `Dm.csv`, but `calc_radar_retrieval_sweep.py` looks for
  `Dm43.csv` (the key in `radarmodule.retrieval_metadata`).

### `DSDlib.retrieval_Cao`
- Uses `radar.K2` computed from the water permittivity, while `polarimetric.calpolrain*` uses
  `Kw2 = 0.93`, a small inconsistency between the retrieval and the forward operator.
- A comment says the lambda search range is in cm^-1; the values are actually mm^-1.
- Fall speeds in the retrieval have no air-density correction.

### Air-density corrections to fall speed are inconsistent across pyPIPS
- `PIPS.calc_empirical_fallspeed_rain`: (1.204/rho)^0.4; `DSDlib.calc_evap`: (1.225/rho)^0.5;
  `DSDlib.calc_VQR_ferrier`: (1.204/rho)^0.5; `PIPS.calc_empirical_fallspeed_hail`:
  (1.225/rho)^0.5. The new `pyPIPS.evaporation` fall-speed classes make rho0 and the exponent
  explicit parameters with defaults matching each existing function; no existing default was
  changed.
