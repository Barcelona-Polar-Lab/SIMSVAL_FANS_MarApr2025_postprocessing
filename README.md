# ARICE 2025 — CTD Post-Processing Pipeline

MATLAB pipeline for post-processing CTD data from **RBR Concerto** sensors deployed during the **ARICE-PONANT 2025** spring campaigns (*SIMSVAL* and *FANS*) aboard R/V *Le Commandant Charcot* in western Greenland fjords (60.2–69.5°N).

This code produces the datasets published in:

> Hoareau et al. (2025). *Hydrographic CTD profiles capturing the onset of near-surface stratification in West Greenland fjords during the SIMSVAL and FANS spring 2025 campaigns.* Earth System Science Data. *(in preparation)*

---

## Dataset overview

| | Value |
|---|---|
| Region | Western Greenland fjords |
| Campaigns | SIMSVAL (Mar–Apr 2025) + FANS (Apr 2025) |
| Sensors | RBR Concerto³ 16 Hz (ICM-CSIC) + RBR Concerto 2 Hz (PONANT) |
| L2 profiles | 73 (32 SIMSVAL + 41 FANS) |
| L0 profiles | 80 |

---

## Repository structure

```
ARICE2025-CTD-postprocessing/
├── scripts/                          ← MATLAB processing scripts
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

## Requirements

- **MATLAB** R2021a or later
- **RSKtools** v3.6 — [RBR Ltd.](https://rbr-global.com/support/matlab-tools/)
- **TEOS-10 GSW** toolbox — [TEOS-10.org](http://www.teos-10.org/software.htm)
- **cmocean** colormaps *(optional, for figures)* — [matplotlib.org/cmocean](https://matplotlib.org/cmocean/)

---

## How to run

### L0 — Raw NetCDF (no processing)

```matlab
% Step 1a: Extract raw SIMSVAL profiles from .rsk files
run('scripts/raw_concat_SIMSVAL_matfile.m')

% Step 1b: Extract raw FANS profiles from .rsk files
run('scripts/raw_concat_FANS_matfile.m')

% Step 2: Combine and export → RAW_ARICE-2025_Greenland_CTD.nc
run('scripts/raw_run_CTD_concat_oceancasts_NetCDF_export.m')
```

### L2 — Processed NetCDF (main product)

```matlab
% Step 1: Process each station individually
% Edit proc_run_CTD_by_stations.m: set campaign_type and station, then run
run('scripts/proc_run_CTD_by_stations.m')

% Step 2a: Concatenate SIMSVAL stations → PROC_CTD_SIMSVAL_oceanCasts.mat
run('scripts/proc_concat_SIMSVAL_matfile.m')

% Step 2b: Concatenate FANS stations → PROC_CTD_FANS_oceanCasts.mat
run('scripts/proc_concat_FANS_matfile.m')

% Step 3: Combine and export → PROC_CTD_ARICE_2025_Greenland_oceanCasts.nc
run('scripts/proc_run_CTD_concat_oceancasts_NetCDF_export.m')
```

> **Note:** Raw `.rsk` data files are not distributed in this repository (too large). They are archived separately and available upon request. Scripts use relative paths — run from the `scripts/` directory or adjust `addpath` calls for your system.

---

## Processing pipeline (L2)

| Step | Operation |
|------|-----------|
| 1 | Sea pressure correction (Patm from FerryBox) |
| 2 | A/D hold correction — 16 Hz only (`RSKcorrecthold`) |
| 3 | Despiking — 4σ, window 15 pts (`RSKdespike`) |
| 4 | C/T lag correction — 16 Hz only, cap ±2 scans (`RSKalignchannel`) |
| 5 | Smoothing — window 5 pts (`RSKsmooth`) |
| 6 | Depth + velocity derivation |
| 7 | Soak removal — fixed_time 20 s (`RSKtrim_soak`) |
| 8 | Loop removal — threshold 0.1 m/s (`RSKremoveloops`) |
| 9 | Salinity + σ-θ derivation (TEOS-10) |
| 10 | QC flags — 7 automated tests, SeaDataNet L20 (`apply_QC_tests`) |

---

## Output files

| File | Level | Profiles | Description |
|------|-------|----------|-------------|
| `RAW_ARICE-2025_Greenland_CTD.nc` | L0 | 80 | Raw data, no processing |
| `PROC_CTD_ARICE_2025_Greenland_oceanCasts.nc` | L2 | 73 | Processed + QC |

Both files follow **CF-1.8 / ACDD-1.3** conventions, `featureType = "profile"`, 2D NaN-padded structure `(obs × profile)`.

---

## Authors

Nina Hoareau, Maria Sánchez, Júlia Crespin, Eva De-Andrés, Ferran Hernández-Macià, Carolina Gabarró, Marta Umbert

Institut de Ciències del Mar (ICM-CSIC), Barcelona, Spain

Contact: nhoareau@icm.csic.es

---

## License

MIT — see [LICENSE](LICENSE)

---

## Citation

If you use this code, please cite:

```bibtex
@software{hoareau2025_ctd_pipeline,
  author    = {Hoareau, Nina and Sánchez, Maria and Crespin, Júlia and
               De-Andrés, Eva and Hernández-Macià, Ferran and
               Gabarró, Carolina and Umbert, Marta},
  title     = {ARICE 2025 CTD post-processing pipeline for RBR Concerto sensors},
  year      = {2025},
  version   = {1.0.0},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.XXXXXXX}
}
```
