%% Create a combined file with all ocean CTD data
% SIMSVAL Greenland - ARICE 2025
% RBR concerto3
% VERSION 5: with QC flags, soak info, validation and qc_summary
% ...
% Nina Hoareau, Oct. 2025
% ...

clear all
close all

%% Load station files
simsval_path = '../outputs/';

Adata = load([simsval_path 'PROC_CTD_SIMSVAL_A.mat']);
Bdata = load([simsval_path 'PROC_CTD_SIMSVAL_B.mat']);
Cdata = load([simsval_path 'PROC_CTD_SIMSVAL_C.mat']);
Idata = load([simsval_path 'PROC_CTD_SIMSVAL_I.mat']);
bd    = load([simsval_path 'PROC_CTD_SIMSVAL_B1_val.mat']);

% Station summary (fixed names, counts from loaded data)
stn_summary = struct( ...
    'name', {'A',                 'B',                 'C',                 'I',                 'B1_val'}, ...
    'n',    {length(Adata.ctdid), length(Bdata.ctdid), length(Cdata.ctdid), length(Idata.ctdid), length(bd.ctdid)});

%% Concatenate atmospheric variables
atmhPa1 = [Adata.atmhPa1 Bdata.atmhPa1 Cdata.atmhPa1 Idata.atmhPa1 bd.atmhPa1]
atmhPa2 = [Adata.atmhPa2 Bdata.atmhPa2 Cdata.atmhPa2 Idata.atmhPa2 bd.atmhPa2]
atmrelativehumidity = [Adata.atmrelativehumidity Bdata.atmrelativehumidity ...
    Cdata.atmrelativehumidity Idata.atmrelativehumidity bd.atmrelativehumidity]
atmwinddirection = [Adata.atmwinddirection Bdata.atmwinddirection ...
    Cdata.atmwinddirection Idata.atmwinddirection bd.atmwinddirection]
atmwindspeed = [Adata.atmwindspeed Bdata.atmwindspeed Cdata.atmwindspeed ...
    Idata.atmwindspeed bd.atmwindspeed]

%% Concatenate metadata
ctdid = [Adata.ctdid Bdata.ctdid Cdata.ctdid Idata.ctdid bd.ctdid]
% ctddate may be row or column vector depending on the file — force row
ctddate = [Adata.ctddate(:)' Bdata.ctddate(:)' Cdata.ctddate(:)' Idata.ctddate(:)' bd.ctddate(:)']
ctdlat = [Adata.ctdlat(:)' Bdata.ctdlat(:)' Cdata.ctdlat(:)' Idata.ctdlat(:)' bd.ctdlat(:)']
ctdlon = [Adata.ctdlon(:)' Bdata.ctdlon(:)' Cdata.ctdlon(:)' Idata.ctdlon(:)' bd.ctdlon(:)']

%% NaN-padding — align all stations to obs_max rows
all_stations = {Adata, Bdata, Cdata, Idata, bd};
obs_max = max(cellfun(@(d) size(d.ctdtemp, 1), all_stations));
fprintf('  obs_max (rows): %d\n', obs_max);

pad    = @(M) [M; NaN(obs_max - size(M,1), size(M,2))];
pad_qc = @(M) [M; 9*ones(obs_max - size(M,1), size(M,2), 'int8')];

%% Concatenate CTD data
ctdseapres = [pad(Adata.ctdseapres) pad(Bdata.ctdseapres) pad(Cdata.ctdseapres) ...
    pad(Idata.ctdseapres) pad(bd.ctdseapres)];

ctdsampltime = [pad(Adata.ctdsampltime) pad(Bdata.ctdsampltime) pad(Cdata.ctdsampltime) ...
    pad(Idata.ctdsampltime) pad(bd.ctdsampltime)];

ctddepth = [pad(Adata.ctddepth) pad(Bdata.ctddepth) pad(Cdata.ctddepth) ...
    pad(Idata.ctddepth) pad(bd.ctddepth)];
ctddepth_mean = mean(ctddepth, 2, 'omitnan');

ctdcond = [pad(Adata.ctdcond) pad(Bdata.ctdcond) pad(Cdata.ctdcond) ...
    pad(Idata.ctdcond) pad(bd.ctdcond)];
ctdtemp = [pad(Adata.ctdtemp) pad(Bdata.ctdtemp) pad(Cdata.ctdtemp) ...
    pad(Idata.ctdtemp) pad(bd.ctdtemp)];
ctdsal = [pad(Adata.ctdsal) pad(Bdata.ctdsal) pad(Cdata.ctdsal) ...
    pad(Idata.ctdsal) pad(bd.ctdsal)];
ctdsigma = [pad(Adata.ctdsigma) pad(Bdata.ctdsigma) pad(Cdata.ctdsigma) ...
    pad(Idata.ctdsigma) pad(bd.ctdsigma)];
ctdvelprof = [pad(Adata.ctdvelprof) pad(Bdata.ctdvelprof) pad(Cdata.ctdvelprof) ...
    pad(Idata.ctdvelprof) pad(bd.ctdvelprof)];

%% Concatenate QC flags
ctdseapres_qc = [pad_qc(Adata.ctdseapres_qc) pad_qc(Bdata.ctdseapres_qc) ...
    pad_qc(Cdata.ctdseapres_qc) pad_qc(Idata.ctdseapres_qc) pad_qc(bd.ctdseapres_qc)];
ctdtemp_qc = [pad_qc(Adata.ctdtemp_qc) pad_qc(Bdata.ctdtemp_qc) ...
    pad_qc(Cdata.ctdtemp_qc) pad_qc(Idata.ctdtemp_qc) pad_qc(bd.ctdtemp_qc)];
ctdcond_qc = [pad_qc(Adata.ctdcond_qc) pad_qc(Bdata.ctdcond_qc) ...
    pad_qc(Cdata.ctdcond_qc) pad_qc(Idata.ctdcond_qc) pad_qc(bd.ctdcond_qc)];
ctdsal_qc = [pad_qc(Adata.ctdsal_qc) pad_qc(Bdata.ctdsal_qc) ...
    pad_qc(Cdata.ctdsal_qc) pad_qc(Idata.ctdsal_qc) pad_qc(bd.ctdsal_qc)];
ctdsigma_qc = [pad_qc(Adata.ctdsigma_qc) pad_qc(Bdata.ctdsigma_qc) ...
    pad_qc(Cdata.ctdsigma_qc) pad_qc(Idata.ctdsigma_qc) pad_qc(bd.ctdsigma_qc)];

%% Concatenate soak info
soak_duration_s = [Adata.soak_duration_s Bdata.soak_duration_s Cdata.soak_duration_s ...
    Idata.soak_duration_s bd.soak_duration_s]
soak_depth_dbar = [Adata.soak_depth_dbar Bdata.soak_depth_dbar Cdata.soak_depth_dbar ...
    Idata.soak_depth_dbar bd.soak_depth_dbar]
soak_n_filtered = [Adata.soak_n_filtered Bdata.soak_n_filtered Cdata.soak_n_filtered ...
    Idata.soak_n_filtered bd.soak_n_filtered]

%% Save
outfile = [simsval_path 'PROC_CTD_SIMSVAL_oceanCasts.mat'];
save(outfile, ...
    'ctddate', 'ctdid', 'ctdlat', 'ctdlon', 'ctddepth_mean', 'ctdseapres', 'ctdsampltime', ...
    'ctddepth', 'ctdcond', 'ctdtemp', 'ctdsal', 'ctdsigma', 'ctdvelprof', ...
    'ctdseapres_qc', 'ctdtemp_qc', 'ctdcond_qc', 'ctdsal_qc', 'ctdsigma_qc', ...
    'soak_duration_s', 'soak_depth_dbar', 'soak_n_filtered', ...
    'atmhPa1', 'atmhPa2', 'atmwindspeed', 'atmwinddirection', 'atmrelativehumidity');

disp(' ')
disp('════════════════════════════════════════════════════════════════')
disp('  CONCATENATION COMPLETE - VERSION 5')
disp('════════════════════════════════════════════════════════════════')
fprintf('  File saved: %s\n', outfile)
fprintf('  Total profiles: %d\n', length(ctdid))
disp(' ')
disp('  Stations included:')
for i = 1:length(stn_summary)
    fprintf('    - %s (%d profile%s)\n', stn_summary(i).name, stn_summary(i).n, ...
        repmat('s', 1, stn_summary(i).n > 1));
end
disp(' ')
disp('  Variables included:')
disp('    - CTD data: ctdtemp, ctdsal, ctdcond, ctdsigma, ctdvelprof')
disp('    - QC flags: ctdtemp_qc, ctdsal_qc, ctdcond_qc, ctdsigma_qc, ctdseapres_qc')
disp('    - Soak info: soak_duration_s, soak_depth_dbar, soak_n_filtered')
disp('    - Atmosphere: atmhPa1, atmhPa2, atmwindspeed, atmwinddirection, atmrelativehumidity')
disp('════════════════════════════════════════════════════════════════')
