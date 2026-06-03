%% concat_raw_SIMSVAL_matfile.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nina HOAREAU, Mai 2026 (BPL, ICM-CSIC)
% ---
% Extraction of raw CTD data (unprocessed) from .rsk files
% for the SIMSVAL campaign (ARICE 2025).
%
% .rsk sources (each station in its own file):
%   A       (16 Hz) : 20250324_A.rsk        — 8 profils
%   B       (16 Hz) : 20250324_B.rsk        — 3 profils
%   B1_val  ( 2 Hz) : 20250324_B1_val.rsk   — 1 profil
%   C       ( 2 Hz) : 20250326_C.rsk        — 11 profils
%   I       (16 Hz) : 237957_20250404_I.rsk — 9 profils
%
% Stations excluded (sea ice): D, G, H
%
% Output : RAW_CTD_SIMSVAL_oceanCasts.mat (cell arrays, 1×32 profils)
%
% IMPORTANT: raw data with no processing whatsoever.
%             Each .rsk file contains ALL downcasts for the station.
%             RSKreadprofiles is called without a specific profile index.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% ════════════════════════════════════════════════════════════════════════
%  PATHS AND TOOLBOXES
%  ════════════════════════════════════════════════════════════════════════

% macOS
addpath('/Users/ninah/Desktop/WORK/MATLAB_toolbox/rbr-rsktools/');
addpath('/Users/ninah/Desktop/WORK/MATLAB_toolbox/mksqlite/');

% Linux (comment/uncomment)
%addpath('/home/nina/Escritorio/WORK/Programs/MATLAB_toolbox/rbr-rsktools/');
%addpath('/home/nina/Escritorio/WORK/Programs/MATLAB_toolbox/mksqlite/');

base_path = '/Users/ninah/Desktop/WORK/ARICE-PONANT/Data_Arctic/2025/';
%base_path = '/home/nina/Escritorio/WORK/BarcelonaPolarLab/Data/in-situ/ARICE-PONANT/Data Arctic/';

out_path = '../outputs/';
atm_file = '../ancillary_data/Atmospheric_Data_SIMSVAL_hourlyMean.mat';

disp(' ')
disp('╔════════════════════════════════════════════════════════════════════╗')
disp('║  RAW EXTRACTION — CTD SIMSVAL Greenland (ARICE 2025)               ║')
disp('╚════════════════════════════════════════════════════════════════════╝')

%% ════════════════════════════════════════════════════════════════════════
%  LOADING ATMOSPHERIC DATA
%  ════════════════════════════════════════════════════════════════════════
disp('Loading atmospheric data...')
atm_data = load(atm_file);
% Structure : atm_data.hourlyMean.DateTime, AtmPress1_hPa, AtmPress2_hPa,
%             WindSpeed_ms, WindDirection, RelativeHumidity

%% ════════════════════════════════════════════════════════════════════════
%  STATION GROUPS DEFINITION - SIMSVAL
%  Each group = one distinct .rsk file
%  All downcasts in the file are read (no specific profile index)
%  ════════════════════════════════════════════════════════════════════════

grp = struct();

grp(1).name   = 'A';
grp(1).lat    = [60.699583, 60.695667, 60.688550, 60.693117, ...
                 60.695383, 60.694983, 60.692967, 60.673933];
grp(1).lon    = [-46.082883, -46.075933, -46.075333, -46.094367, ...
                 -46.106333, -46.160033, -46.165383, -46.156833];
grp(1).ID     = {'A1','A2','A3','A4','A5','A6','A7','A8'};
grp(1).is_2hz = false;
grp(1).rsk_file = '../RAW/raw_SIMSVAL/20250324_A.rsk';

grp(2).name   = 'B';
grp(2).lat    = [60.779000, 60.786750, 60.796450];
grp(2).lon    = [-45.629200, -45.634583, -45.639683];
grp(2).ID     = {'B1','B2','B3'};
grp(2).is_2hz = false;
grp(2).rsk_file = '../RAW/raw_SIMSVAL//20250324_B.rsk';

grp(3).name   = 'B1_val';
grp(3).lat    = 60.779000;
grp(3).lon    = -45.629200;
grp(3).ID     = {'B1_val'};
grp(3).is_2hz = true;
grp(3).rsk_file = '../RAW/raw_SIMSVAL/20250324_B1_val.rsk';

grp(4).name   = 'C';
grp(4).lat    = [60.294683, 60.298683, 60.302633, 60.346950, 60.346267, ...
                 60.345617, 60.265150, 60.241700, 60.229583, 60.228817, 60.224400];
grp(4).lon    = [-44.004317, -44.006933, -44.008333, -44.064733, -44.078833, ...
                 -44.093250, -44.094050, -44.131533, -44.163150, -44.157800, -44.145083];
grp(4).ID     = {'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11'};
grp(4).is_2hz = true;
grp(4).rsk_file = '../RAW/raw_SIMSVAL/20250326_C.rsk';

grp(5).name   = 'I';
grp(5).lat    = [65.985000, 65.986517, 65.992333, 65.984300, ...
                 65.948283, 65.906017, 65.912317, 65.920100, 65.943533];
grp(5).lon    = [-52.523333, -52.530700, -52.537200, -52.590533, ...
                 -52.644367, -52.816067, -52.822633, -52.825783, -52.722583];
grp(5).ID     = {'I1','I2','I3','I4','I5','I6','I7','I8','I9'};
grp(5).is_2hz = false;
grp(5).rsk_file = '../RAW/raw_SIMSVAL/237957_20250404_I.rsk';

% Fjord names per group
grp(1).fjord = 'Qaqortoq Fjord';    % A
grp(2).fjord = 'Eqaluit Fjord';     % B
grp(3).fjord = 'Eqaluit Fjord';     % B1_val
grp(4).fjord = 'Aappilattoq Fjord'; % C
grp(5).fjord = 'Eternity Fjord';    % I

%% ════════════════════════════════════════════════════════════════════════
%  RAW DATA EXTRACTION
%  ════════════════════════════════════════════════════════════════════════

n_total = sum(arrayfun(@(g) length(g.ID), grp));
raw_T   = cell(1, n_total);
raw_C   = cell(1, n_total);
raw_P   = cell(1, n_total);
raw_t   = cell(1, n_total);
n_samp  = zeros(1, n_total);
freq    = zeros(1, n_total);
ctdid   = cell(1, n_total);
ctdfjord = cell(1, n_total);
ctdlat  = NaN(1, n_total);
ctdlon  = NaN(1, n_total);
ctddate = NaT(1, n_total, 'TimeZone', 'UTC');
atmhPa1          = NaN(1, n_total);
atmhPa2          = NaN(1, n_total);
atmwindspeed     = NaN(1, n_total);
atmwinddirection = NaN(1, n_total);
atmrelativehumidity = NaN(1, n_total);

% Convert atm timestamps to datenum once (avoids timezone conflict)
atm_dn   = datenum(atm_data.hourlyMean.DateTime);
epoch_dn = datenum('1970-01-01 00:00:00');

ip = 0;

for ig = 1:length(grp)
    g = grp(ig);
    fprintf('\n--- Station %s ---\n', g.name)
    fprintf('  Opening: %s\n', g.rsk_file)

    rsk = RSKopen(g.rsk_file);

    % Read all downcasts in the file (no specific index — 1 file = 1 station)
    rsk_g = RSKreadprofiles(rsk, 'direction', 'down');
    np_g  = length(rsk_g.data);
    fprintf('  %d profiles read\n', np_g)

    % Consistency check
    n_expected = length(g.ID);
    if np_g ~= n_expected
        warning('Station %s: %d profiles read, %d expected (defined IDs)', g.name, np_g, n_expected)
    end

    % Channel indices
    if g.is_2hz
        try
            idx_T = getchannelindex(rsk_g, 'Temperature1');
        catch
            idx_T = getchannelindex(rsk_g, 'Temperature');
            warning('Station %s: Temperature1 not found, using Temperature', g.name);
        end
        fs = 2;
    else
        idx_T = getchannelindex(rsk_g, 'Temperature');
        fs = 16;
    end
    idx_C = getchannelindex(rsk_g, 'Conductivity');
    idx_P = getchannelindex(rsk_g, 'Pressure');

    for k = 1:np_g
        ip = ip + 1;

        T_k = rsk_g.data(k).values(:, idx_T);
        C_k = rsk_g.data(k).values(:, idx_C);
        P_k = rsk_g.data(k).values(:, idx_P);
        t_k = rsk_g.data(k).tstamp;

        % Timestamps in POSIX seconds (datenum → s since 1970-01-01)
        t_posix_k = (t_k - epoch_dn) * 86400;

        raw_T{ip} = T_k;
        raw_C{ip} = C_k;
        raw_P{ip} = P_k;
        raw_t{ip} = t_posix_k;
        n_samp(ip) = length(T_k);
        freq(ip)   = fs;

        ctdid{ip}   = g.ID{k};
        ctdfjord{ip} = g.fjord;
        ctdlat(ip)  = g.lat(k);
        ctdlon(ip)  = g.lon(k);

        % Profile start time → UTC datetime for ctddate
        ctddate(ip) = datetime(t_k(1), 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');

        % Atmospheric data: match by datenum (avoids timezone conflict)
        [~, iatm] = min(abs(atm_dn - t_k(1)));
        atmhPa1(ip)          = atm_data.hourlyMean.AtmPress1_hPa(iatm);
        atmhPa2(ip)          = atm_data.hourlyMean.AtmPress2_hPa(iatm);
        atmwindspeed(ip)     = atm_data.hourlyMean.WindSpeed_ms(iatm);
        atmwinddirection(ip) = atm_data.hourlyMean.WindDirection(iatm);
        atmrelativehumidity(ip) = atm_data.hourlyMean.RelativeHumidity(iatm);

        fprintf('  [%2d/%2d] %s : %d pts, Patm=%.1f hPa\n', ...
            ip, n_total, g.ID{k}, n_samp(ip), atmhPa1(ip))
    end

    clear rsk rsk_g
end

fprintf('\n  Total: %d profiles extracted\n', ip)
if ip ~= n_total
    warning('Expected %d profiles, got %d', n_total, ip)
end

%% ════════════════════════════════════════════════════════════════════════
%  SAVE
%  ════════════════════════════════════════════════════════════════════════

outfile = [out_path 'RAW_CTD_SIMSVAL_oceanCasts.mat'];
fprintf('\nSaving to: %s\n', outfile)

save(outfile, ...
    'raw_T', 'raw_C', 'raw_P', 'raw_t', ...
    'n_samp', 'freq', ...
    'ctdid', 'ctdfjord', 'ctdlat', 'ctdlon', 'ctddate', ...
    'atmhPa1', 'atmhPa2', 'atmwindspeed', 'atmwinddirection', 'atmrelativehumidity');

d = dir(outfile);
fprintf('  [OK] Saved: %.1f KB\n', d.bytes/1024)

disp(' ')
disp('╔════════════════════════════════════════════════════════════════════╗')
disp('║  SIMSVAL RAW — DONE                                                ║')
disp('╚════════════════════════════════════════════════════════════════════╝')
summary_parts = arrayfun(@(g) sprintf('%s(%d)', g.name, length(g.ID)), grp, 'UniformOutput', false);
fprintf('  %d profils | %s\n', n_total, strjoin(summary_parts, ' + '))
fprintf('  n_samp min/max: %d / %d\n', min(n_samp), max(n_samp))
fprintf('  Frequencies: 16 Hz (A, B, I) | 2 Hz (B1_val, C)\n')
