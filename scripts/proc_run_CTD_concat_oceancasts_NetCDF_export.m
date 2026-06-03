%% Main script: Concatenation and NetCDF export of CTD data
% ARICE 2025 - Greenland Ocean Casts
% Nina Hoareau, Jan. 2025
%
% This script orchestrates:
%   1. Concatenation of CTD data (SIMSVAL + FANS)
%   2. Saving the combined .mat file
%   3. Export to NetCDF (CF-1.8 / ACDD-1.3)
%
% PREREQUISITES:
%   - Per-campaign files must exist in ../outputs/:
%       PROC_CTD_SIMSVAL_oceanCasts.mat & PROC_CTD_FANS_oceanCasts.mat
%   - TEOS-10 GSW toolbox in path
%   - export_netcdf.m in the same folder

clear all
close all

% Toolboxes
% linux path
addpath('/home/nina/Escritorio/WORK/Programs/MATLAB_toolbox/rbr-rsktools/');
addpath('/home/nina/Escritorio/WORK/Programs/MATLAB_toolbox/teos10/');
addpath('/home/nina/Escritorio/WORK/Programs/MATLAB_toolbox/teos10/library/');
addpath('/home/nina/Escritorio/WORK/Programs/MATLAB_toolbox/cmocean/');
% macos path
% addpath('/Users/ninah/Desktop/WORK/MATLAB_toolbox/teos10/');
% addpath('/Users/ninah/Desktop/WORK/MATLAB_toolbox/teos10/library/');
% addpath('/Users/ninah/Desktop/WORK/MATLAB_toolbox/teos10/thermodynamics_from_t/');
% addpath('/Users/ninah/Desktop/WORK/MATLAB_toolbox/cmocean/');

% Output files
matfile = '../outputs/PROC_CTD_ARICE_2025_Greenland_oceanCasts.mat';
ncfile  = '../outputs/PROC_CTD_ARICE_2025_Greenland_oceanCasts.nc';

%% ════════════════════════════════════════════════════════════════════════
%  STEP 1: Concatenation of CTD data (SIMSVAL + FANS)
%  ════════════════════════════════════════════════════════════════════════
disp(' ')
disp('╔════════════════════════════════════════════════════════════════╗')
disp('║  STEP 1: CONCATENATION OF CTD DATA                             ║')
disp('╚════════════════════════════════════════════════════════════════╝')

% Load per-campaign files
disp('Loading campaign files...')
A = load('../outputs/PROC_CTD_SIMSVAL_oceanCasts.mat');
B = load('../outputs/PROC_CTD_FANS_oceanCasts.mat');

fprintf('  SIMSVAL : %d profiles loaded\n', length(A.ctdid))
fprintf('  FANS    : %d profiles loaded\n', length(B.ctdid))

% Concatenate atmospheric variables
atmhPa1 = [A.atmhPa1 B.atmhPa1];
atmhPa2 = [A.atmhPa2 B.atmhPa2];
atmrelativehumidity = [A.atmrelativehumidity B.atmrelativehumidity];
atmwinddirection = [A.atmwinddirection B.atmwinddirection];
atmwindspeed = [A.atmwindspeed B.atmwindspeed];

% Concatenate metadata
ctdid = [A.ctdid B.ctdid];
ctddate = [A.ctddate(:)' B.ctddate(:)'];
ctdlat = [A.ctdlat(:)' B.ctdlat(:)'];
ctdlon = [A.ctdlon(:)' B.ctdlon(:)'];

% Campaign labels
campaign = [repmat({'SIMSVAL'}, 1, length(A.ctdid)), repmat({'FANS'}, 1, length(B.ctdid))];

% NaN-padding — align SIMSVAL and FANS to a common obs_max
% (required at native resolution: obs_max differs between campaigns)
obs_max = max(size(A.ctdtemp, 1), size(B.ctdtemp, 1));
fprintf('  common obs_max SIMSVAL+FANS: %d rows\n', obs_max);

pad    = @(M) [M; NaN(obs_max - size(M,1), size(M,2))];
pad_qc = @(M) [M; 9*ones(obs_max - size(M,1), size(M,2), 'int8')];

% ctdseapres, ctdsampltime and ctddepth — 2D NaN-padded
ctdseapres   = [pad(A.ctdseapres)   pad(B.ctdseapres)];
ctdsampltime = [pad(A.ctdsampltime) pad(B.ctdsampltime)];
mean_lat   = mean(ctdlat, 'omitnan');
ctddepth   = abs(gsw_z_from_p(ctdseapres, mean_lat));
fprintf('  Depth computed at mean latitude = %.2f°N\n', mean_lat)

% Concatenate CTD data
ctdcond    = [pad(A.ctdcond)    pad(B.ctdcond)];
ctdtemp    = [pad(A.ctdtemp)    pad(B.ctdtemp)];
ctdsal     = [pad(A.ctdsal)     pad(B.ctdsal)];
ctdsigma   = [pad(A.ctdsigma)   pad(B.ctdsigma)];
ctdvelprof = [pad(A.ctdvelprof) pad(B.ctdvelprof)];

% Concatenate QC flags
ctdseapres_qc = [pad_qc(A.ctdseapres_qc) pad_qc(B.ctdseapres_qc)];
ctdtemp_qc    = [pad_qc(A.ctdtemp_qc)    pad_qc(B.ctdtemp_qc)];
ctdcond_qc    = [pad_qc(A.ctdcond_qc)    pad_qc(B.ctdcond_qc)];
ctdsal_qc     = [pad_qc(A.ctdsal_qc)     pad_qc(B.ctdsal_qc)];
ctdsigma_qc   = [pad_qc(A.ctdsigma_qc)   pad_qc(B.ctdsigma_qc)];

% Concatenate soak info
soak_duration_s = [A.soak_duration_s B.soak_duration_s];
soak_depth_dbar = [A.soak_depth_dbar B.soak_depth_dbar];
soak_n_filtered = [A.soak_n_filtered B.soak_n_filtered];

fprintf('  Total: %d profiles concatenated\n', length(ctdid))
disp('  [OK] Concatenation done')

%% ════════════════════════════════════════════════════════════════════════
%  STEP 2: Save combined .mat file
%  ════════════════════════════════════════════════════════════════════════
disp(' ')
disp('╔════════════════════════════════════════════════════════════════╗')
disp('║  STEP 2: SAVE .MAT FILE                                        ║')
disp('╚════════════════════════════════════════════════════════════════╝')

save(matfile, ...
    'ctddate', 'ctdid', 'ctdlat', 'ctdlon', 'ctdseapres', 'ctdsampltime', 'ctddepth', 'mean_lat', ...
    'ctdcond', 'ctdtemp', 'ctdsal', 'ctdsigma', 'ctdvelprof', ...
    'ctdseapres_qc', 'ctdtemp_qc', 'ctdcond_qc', 'ctdsal_qc', 'ctdsigma_qc', ...
    'soak_duration_s', 'soak_depth_dbar', 'soak_n_filtered', ...
    'campaign', ...
    'atmhPa1', 'atmhPa2', 'atmwindspeed', 'atmwinddirection', 'atmrelativehumidity');

d = dir(matfile);
fprintf('  File saved: %s\n', matfile)
fprintf('  Size: %.1f KB\n', d.bytes/1024)
disp('  [OK] .mat file saved')

%% ════════════════════════════════════════════════════════════════════════
%  STEP 3: Export NetCDF (CF-1.8 / ACDD-1.3)
%  ════════════════════════════════════════════════════════════════════════
disp(' ')
disp('╔════════════════════════════════════════════════════════════════╗')
disp('║  STEP 3: NETCDF EXPORT                                         ║')
disp('╚════════════════════════════════════════════════════════════════╝')

export_netcdf(matfile, ncfile, 'sss_file', '../ancillary_data/SSS_autosal.csv');

%% ════════════════════════════════════════════════════════════════════════
%  FINAL SUMMARY
%  ════════════════════════════════════════════════════════════════════════
disp(' ')
disp('╔════════════════════════════════════════════════════════════════╗')
disp('║  SUMMARY - CTD Greenland - ARICE 2025                          ║')
disp('╚════════════════════════════════════════════════════════════════╝')
disp(' ')
disp('  Generated files:')
fprintf('    1. %s (.mat v4)\n', matfile)
fprintf('    2. %s (NetCDF-4 CF-1.8)\n', ncfile)
disp(' ')
fprintf('  Number of profiles: %d\n', length(ctdid))
fprintf('  Campaigns: SIMSVAL (%d) + FANS (%d)\n', length(A.ctdid), length(B.ctdid))
disp(' ')
disp('  CTD variables:')
disp('    - Temperature, Salinity, Conductivity, Sigma-theta, Velocity')
disp('    - Soak info (duration, depth, filtered points)')
disp('    - Atmospheric data (pressure, wind, humidity)')
disp(' ')
disp('  Quality Control:')
disp('    - Type: Automatic (RTQC-like, QARTOD-style)')
disp('    - 7 tests applied:')
disp('        1. NaN check')
disp('        2. Gross range (global)')
disp('        3. Regional range (Arctic)')
disp('        4. Flat line')
disp('        5. Vertical gradient')
disp('        6. Density inversion')
disp('        7. Pressure monotonicity')
disp('    - Flags SeaDataNet : 1=good, 2=probably good, 3=suspect, 4=bad, 9=missing')
disp('    - NOTE: No DMQC validation (no comparison with climatologies/floats)')
disp(' ')
disp('  Processing level:')
disp('    - L2: Post-processed with automated QC (no delayed-mode validation)')
disp(' ')
disp('  Next steps:')
disp('    - CF validation: cfchecks PROC_CTD_ARICE_2025_Greenland_oceanCasts.nc')
disp('    - Inspection: ncdisp(''PROC_CTD_ARICE_2025_Greenland_oceanCasts.nc'')')
disp('════════════════════════════════════════════════════════════════════')
