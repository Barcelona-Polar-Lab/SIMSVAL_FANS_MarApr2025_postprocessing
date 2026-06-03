# SIMSVAL & FANS Mar–Apr 2025 — CTD Post-Processing Pipeline

[![DOI code](https://img.shields.io/badge/DOI%20code-10.5281%2Fzenodo.20527364-blue)](https://doi.org/10.5281/zenodo.20527364)
[![DOI data](https://img.shields.io/badge/DOI%20data-10.5281%2Fzenodo.20489479-blue)](https://doi.org/10.5281/zenodo.20489479)

MATLAB pipeline for post-processing CTD data from **RBR Concerto** sensors deployed during the **SIMSVAL and FANS** 2025 spring campaigns aboard R/V *Le Commandant Charcot* in western Greenland fjords (60.2–69.5°N).

This code produces the datasets published in:

> Umbert et al. (2026). *Hydrographic CTD profiles capturing the onset of near-surface stratification in Southwest Greenland fjords during two consecutive spring 2025 field campaigns.* Earth System Science Data. *(in preparation)*

---

## Dataset overview

| | Value |
|---|---|
| Region | Western Greenland fjords |
| Campaigns | SIMSVAL (Mar–Apr 2025) + FANS (Apr 2025) |
| Sensors | RBR Concerto³ 16 Hz (ICM-CSIC, S/N 237957) + RBR Concerto 2 Hz (PONANT) |
| L2 profiles | 73 (32 SIMSVAL + 41 FANS) |
| L0 profiles | 73 |

---

## Requirements

### MATLAB version
MATLAB R2021a or later.

### Required toolboxes

| Toolbox | Version tested | Download |
|---------|---------------|----------|
| **RSKtools** (RBR Ltd.) | v3.6 | [github.com/RBRglobal/rbr-rsktools](https://github.com/RBRglobal/rbr-rsktools) |
| **TEOS-10 GSW** | v3.06.16 | [teos-10.org/software.htm](http://www.teos-10.org/software.htm) |
| **cmocean** *(optional, figures only)* | v2.0 | [github.com/chadagreene/cmocean](https://github.com/chadagreene/cmocean) |

> **Note:** RSKtools requires `mksqlite` to read `.rsk` files. A compiled binary is included with RSKtools for Windows, macOS, and Linux.

### Path configuration

After downloading the toolboxes, run `setup.m` to configure your paths interactively:

```matlab
run('setup.m')
```

Or add them manually to your MATLAB path (edit paths to match your system):

```matlab
addpath('/path/to/rbr-rsktools/');
addpath('/path/to/teos10/');
addpath('/path/to/teos10/library/');
addpath('/path/to/cmocean/');   % optional
```

---

## Raw data

The raw `.rsk` files (RBR binary format) are **not distributed in this repository** due to file size.  
They are available on Zenodo ([doi:10.5281/zenodo.20489479](https://doi.org/10.5281/zenodo.20489479), restricted access — contact nhoareau@icm.csic.es):

| Archive | Contents | Size |
|---------|----------|------|
| `raw_SIMSVAL_rsk.zip` | SIMSVAL .rsk files (stations A, B, B1_val, C, I) | ~8 MB |
| `raw_FANS_rsk.zip` | FANS .rsk files (ICM-CSIC 16 Hz, S/N 237957) | ~16 MB |

Download the `.rsk` files and place them in:
```
scripts/RAW/raw_SIMSVAL/    ← SIMSVAL .rsk files (stations A, B, B1_val, C, I)
scripts/RAW/raw_FANS/       ← FANS .rsk files (ICM 16 Hz)
```

---

## Repository structure

```
SIMSVAL_FANS_MarApr2025_postprocessing/
├── setup.m                           ← Configure MATLAB toolbox paths
├── README.md
├── LICENSE
├── CITATION.cff
├── scripts/
│   ├── proc_run_CTD_by_stations.m    ← Entry point: L2 per-station processing
│   ├── process_RBR_CTD.m             ← Core processing function
│   ├── proc_concat_SIMSVAL_matfile.m ← Concatenate SIMSVAL stations → .mat
│   ├── proc_concat_FANS_matfile.m    ← Concatenate FANS stations → .mat
│   ├── proc_run_CTD_concat_oceancasts_NetCDF_export.m  ← Final L2 NetCDF export
│   ├── raw_concat_SIMSVAL_matfile.m  ← Extract RAW SIMSVAL → .mat
│   ├── raw_concat_FANS_matfile.m     ← Extract RAW FANS → .mat
│   ├── raw_run_CTD_concat_oceancasts_NetCDF_export.m   ← Final L0 NetCDF export
│   ├── RSKtrim_soak.m                ← Soak detection (fixed_time / velocity)
│   ├── apply_QC_tests.m              ← 7 QC tests (SeaDataNet L20 flags)
│   └── export_netcdf.m               ← NetCDF-4 export (CF-1.8 / ACDD-1.3)
└── ancillary_data/
    ├── SSS_autosal.csv               ← AutoSal surface salinity water samples
    ├── station_coordinates.csv       ← Station lat/lon/campaign/frequency/fjord
    ├── Atmospheric_Data_SIMSVAL_hourlyMean.mat  ← FerryBox atmospheric data
    └── Atmospheric_Data_FANS_hourlyMean.mat     ← FerryBox atmospheric data
```

---

## How to run

All scripts use **relative paths** — run from the `scripts/` directory in MATLAB.

### L0 — Raw NetCDF (no processing)

```matlab
cd scripts/

% Step 1a — Extract raw SIMSVAL profiles from .rsk files
run('raw_concat_SIMSVAL_matfile.m')

% Step 1b — Extract raw FANS profiles from .rsk files
run('raw_concat_FANS_matfile.m')

% Step 2 — Combine and export → SIMSVAL_FANS_MarApr2025_RAW.nc
run('raw_run_CTD_concat_oceancasts_NetCDF_export.m')
```

### L2 — Processed NetCDF (main product)

```matlab
cd scripts/

% Step 1 — Process each station individually
% Edit proc_run_CTD_by_stations.m: set campaign_type ('SIMSVAL' or 'FANS')
% and station ('A', 'B', ..., '1A', '1B', ...) then run
run('proc_run_CTD_by_stations.m')

% Step 2a — Concatenate SIMSVAL → PROC_CTD_SIMSVAL_oceanCasts.mat
run('proc_concat_SIMSVAL_matfile.m')

% Step 2b — Concatenate FANS → PROC_CTD_FANS_oceanCasts.mat
run('proc_concat_FANS_matfile.m')

% Step 3 — Combine and export → SIMSVAL_FANS_MarApr2025_PROC.nc
run('proc_run_CTD_concat_oceancasts_NetCDF_export.m')
```

---

## Processing pipeline (L2)

| Step | Operation | Notes |
|------|-----------|-------|
| 1 | Data reading | `RSKopen` + `RSKreadprofiles` (downcast only). Channel `Temperature1` renamed to `Temperature` for the 2 Hz PONANT sensor (station C). |
| 2 | Sea pressure = Pressure − Patm | Patm from FerryBox hourly means, nearest-neighbour interpolation |
| 3 | A/D hold correction (`RSKcorrecthold`) | 16 Hz only |
| 4 | Despiking (`RSKdespike`) | 4σ, window 15 pts, direction *down*, action `nan` |
| 5 | C/T lag correction (`RSKcalculateCTlag` + `RSKalignchannel`) | 16 Hz only, `seapressurerange = [5 30]`, cap ±2 scans |
| 6 | Smoothing (`RSKsmooth`) | window 5 pts, applied to temperature and conductivity |
| 7 | Depth + velocity derivation | `RSKderivedepth` + `RSKderivevelocity` |
| 8 | Soak removal (`RSKtrim_soak`) | fixed_time 20 s (5 s for 1A_1 and 7H_4) |
| 9 | Loop removal (`RSKremoveloops`) | threshold 0.1 m/s |
| 10 | Salinity + σ-θ derivation | `RSKderivesalinity` + `gsw_sigma0` (TEOS-10 GSW) |
| 11 | QC flags (`apply_QC_tests`) | 7 tests, SeaDataNet L20 flags |

> **Vertical-gradient test note:** The QARTOD/Argo vertical-gradient test (test 5) is **not applied** to native-resolution data. At a typical descent velocity of ~0.3 m s⁻¹, the pressure interval between consecutive 16 Hz samples is ~0.003 dbar, which amplifies gradient noise ~80× compared with 0.25 dbar binned data. This test is reserved for bin-averaged products.

---

## Output files

| File | Level | Profiles | Description |
|------|-------|----------|-------------|
| `SIMSVAL_FANS_MarApr2025_RAW.nc` | L0 | 73 | Raw data, atmospheric pressure correction only |
| `SIMSVAL_FANS_MarApr2025_PROC.nc` | L2 | 73 | Processed + QC |

Both follow **CF-1.8 / ACDD-1.3** conventions, `featureType = "profile"`, 2D NaN-padded `(obs × profile)`.

### Dimensions (`SIMSVAL_FANS_MarApr2025_PROC.nc`)

| Dimension | Size | Description |
|-----------|-----:|-------------|
| `profile` | 73 | One element per CTD profile |
| `obs` | 6623 | Maximum number of samples per profile (NaN-padded); native resolution |
| `name_strlen` | 6 | Maximum length of station identifiers |
| `campaign_strlen` | 7 | Maximum length of campaign names (`SIMSVAL` = 7 characters) |

### Key variables

**Coordinate variables** — per profile: `time`, `latitude`, `longitude`, `station_id`, `campaign`

**Coordinate variables** — per sample `(obs, profile)`: `sea_pressure`, `depth`, `sample_time`

**Oceanographic variables** — `(obs, profile)`: `temperature`, `salinity`, `conductivity`, `sigma_theta`, `profile_velocity`

**QC flags** — `(obs, profile)`, int8, SeaDataNet L20: `temperature_qc`, `salinity_qc`, `conductivity_qc`, `sigma_theta_qc`, `sea_pressure_qc`

**Auxiliary per-profile variables**: `n_samples`, `max_sea_pressure`, `sampling_frequency`, `sss` (AutoSal surface salinity, complementary reference only), soak metadata (`soak_duration`, `soak_depth`, `soak_n_filtered`), FerryBox atmospheric variables

---

## Authors

Nina Hoareau, Marta Umbert, Júlia Crespin, María Sánchez Urrea, Ferran Hernández-Macià, Carolina Gabarró, Eva De Andrés

Institut de Ciències del Mar (ICM-CSIC), Barcelona, Spain — nhoareau@icm.csic.es

**Groups:** [Barcelona Polar Lab](https://www.barcelonapolarlab.es) · Servicio de Ingeniería Oceanográfica — ICM-CSIC

---

## Acknowledgements

This work was supported by the Spanish Ministry of Science through the "Severo Ochoa Centre of Excellence" accreditation awarded to ICM-CSIC (CEX2024-001494-S).

---

## License

MIT — see [LICENSE](LICENSE)

---

## Citation

If you use this code, please cite:

```bibtex
@software{umbert2026_ctd_pipeline,
  author    = {Hoareau, Nina and Umbert, Marta and Crespin, Júlia and
               Sánchez Urrea, María and Hernández-Macià, Ferran and
               Gabarró, Carolina and De Andrés, Eva},
  title     = {CTD post-processing pipeline for RBR Concerto sensors —
               SIMSVAL and FANS Mar-Apr 2025},
  year      = {2026},
  version   = {1.0.0},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.20527364}
}
```
