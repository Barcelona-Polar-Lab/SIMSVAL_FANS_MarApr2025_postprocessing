%% concat_raw_FANS_matfile.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nina HOAREAU, Avril 2026 (BPL, ICM-CSIC)
% ---
% Extraction of raw CTD data (unprocessed) from .rsk files
% for the FANS campaign (ARICE 2025).
%
% .rsk sources:
%   ICM (16 Hz) : 237957_20250429_1235.rsk  — ONE file, profiles indexed
%   JP  ( 2 Hz) : 237329_20250422_1459_JPctd.rsk  — fichier séparé
%
% Station to RSK profile index mapping (cf. process_FANS_station_v3.m):
%   profile 1  = A1 failed cast (excluded)
%   F1A  = profiles 2:6    F1B  = profiles 7:9    F2C  = profiles 10:12
%   ICMctd_20  = profiles 13:14  (intercalated between 2C and 3D)
%   F3D  = profiles 15:19  (profile 20 = transit cast, excluded)
%   F4E  = profiles 21:25  F5F  = profiles 26:30  F6G  = profiles 31:34
%   F7H  = profiles 35:38  F7I  = profiles 39:41  F8J  = profiles 42:46
%   JPctd_19 = profile 1 (JP file)
%   JPctd_20 = profiles 2:4 (JP file)
%
% Output : RAW_CTD_FANS_oceanCasts.mat (cell arrays, 1×48 profils)
%
% IMPORTANT: raw data with no processing whatsoever.
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

out_path  = '../outputs/';
fans_rsk_ICM = '../RAW/raw_FANS/237957_20250429_1235.rsk';
fans_rsk_JP  = '../RAW/raw_FANS/237329_20250422_1459_JPctd.rsk';
atm_file     = '../ancillary_data/Atmospheric_Data_FANS_hourlyMean.mat';

disp(' ')
disp('╔════════════════════════════════════════════════════════════════════╗')
disp('║  RAW EXTRACTION — CTD FANS Greenland (ARICE 2025)                  ║')
disp('╚════════════════════════════════════════════════════════════════════╝')

%% ════════════════════════════════════════════════════════════════════════
%  LOADING ATMOSPHERIC DATA
%  ════════════════════════════════════════════════════════════════════════
disp('Loading atmospheric data...')
atm_data = load(atm_file);
% Structure : atm_data.hourlyMean.DateTime, AtmPress1_hPa, AtmPress2_hPa,
%             WindSpeed_ms, WindDirection, RelativeHumidity

%% ════════════════════════════════════════════════════════════════════════
%  STATION GROUPS DEFINITION - FANS
%  Each group = a set of consecutive profiles within the ICM file
%  Same order as the processed file: F1A→F1B→...→F8J→JP19→JP20→ICM20
%  ════════════════════════════════════════════════════════════════════════

% ICM 16 Hz groups (in fans_rsk_ICM file)
grp = struct();
grp(1).name     = 'F1A';
grp(1).profiles = 2:6;   % profile 1 = A1 failed cast (v3 ref)
grp(1).lat      = [66.9417583, 66.9428383, 66.9447300, 66.9447000, 66.9437967];
grp(1).lon      = [-53.7110900, -53.7285250, -53.7442583, -53.7461067, -53.7461067];
grp(1).ID       = {'1A_1','1A_2','1A_3','1A_4','1A_5'};
grp(1).is_2hz   = false;
grp(1).rsk_file = fans_rsk_ICM;

grp(2).name     = 'F1B';
grp(2).profiles = 7:9;
grp(2).lat      = [66.9419883, 66.9417800, 66.9574667];
grp(2).lon      = [-53.7690100, -53.7628567, -53.7525000];
grp(2).ID       = {'1B_1','1B_2','1B_3'};
grp(2).is_2hz   = false;
grp(2).rsk_file = fans_rsk_ICM;

grp(3).name     = 'F2C';
grp(3).profiles = 10:12;
grp(3).lat      = [66.8891983, 66.8931433, 66.8789650];
grp(3).lon      = [-53.8549033, -53.6419133, -53.6877367];
grp(3).ID       = {'2C_1','2C_2','2C_3'};
grp(3).is_2hz   = false;
grp(3).rsk_file = fans_rsk_ICM;

grp(4).name     = 'F3D';
grp(4).profiles = 15:19;  % profiles 13-14 = ICMctd_20 intercalated between 2C and 3D (v3 ref)
grp(4).lat      = [69.4341750, NaN, 69.4532033, 69.4938617, 69.4328983];
grp(4).lon      = [-50.8942933, NaN, -51.0013517, -51.0517317, -51.0828333];
grp(4).ID       = {'3D_1','3D_2','3D_3','3D_4','3D_5'};
grp(4).is_2hz   = false;
grp(4).rsk_file = fans_rsk_ICM;

grp(5).name     = 'F4E';
grp(5).profiles = 21:25;  % profile 20 = transit cast between 3D and 4E (v3 ref)
grp(5).lat      = [61.1372833, 61.1269733, 61.1270667, 61.1283317, 61.1291317];
grp(5).lon      = [-45.4901250, -45.4807433, -45.4651233, -45.4426867, -45.4210700];
grp(5).ID       = {'4E_1','4E_2','4E_3','4E_4','4E_5'};
grp(5).is_2hz   = false;
grp(5).rsk_file = fans_rsk_ICM;

grp(6).name     = 'F5F';
grp(6).profiles = 26:30;
grp(6).lat      = [61.0162367, 61.0229967, 61.0316933, 61.0431567, 61.0040883];
grp(6).lon      = [-46.1302633, -46.1286933, -46.1269300, -46.1303733, -46.0789850];
grp(6).ID       = {'5F_1','5F_2','5F_3','5F_4','5F_5'};
grp(6).is_2hz   = false;
grp(6).rsk_file = fans_rsk_ICM;

grp(7).name     = 'F6G';
grp(7).profiles = 31:34;
grp(7).lat      = [60.7106683, 60.7056267, 60.7000133, 60.6918383];
grp(7).lon      = [-46.0123300, -46.0086967, -46.0044550, -46.0022920];
grp(7).ID       = {'6G_1','6G_2','6G_3','6G_4'};
grp(7).is_2hz   = false;
grp(7).rsk_file = fans_rsk_ICM;

grp(8).name     = 'F7H';
grp(8).profiles = 35:38;
grp(8).lat      = [60.2873383, 60.2913917, 60.2945233, 60.2967983];
grp(8).lon      = [-44.2383067, -44.2293433, -44.2200617, -44.2125383];
grp(8).ID       = {'7H_1','7H_2','7H_3','7H_4'};
grp(8).is_2hz   = false;
grp(8).rsk_file = fans_rsk_ICM;

grp(9).name     = 'F7I';
grp(9).profiles = 39:41;
grp(9).lat      = [60.2638117, 60.2625383, 60.2603300];
grp(9).lon      = [-44.1738883, -44.1916250, -44.2118333];
grp(9).ID       = {'7I_1','7I_2','7I_3'};
grp(9).is_2hz   = false;
grp(9).rsk_file = fans_rsk_ICM;

grp(10).name     = 'F8J';
grp(10).profiles = 42:46;
grp(10).lat      = [60.1801050, 60.1791633, 60.1782383, 60.1801383, 60.1609800];
grp(10).lon      = [-43.6296067, -43.6093467, -43.6203117, -43.6380100, -43.6284800];
grp(10).ID       = {'8J_1','8J_2','8J_3','8J_4','8J_5'};
grp(10).is_2hz   = false;
grp(10).rsk_file = fans_rsk_ICM;

% % JP stations (2 Hz, separate file)
% grp(11).name     = 'JPctd_19';
% grp(11).profiles = 1;
% grp(11).lat      = NaN;
% grp(11).lon      = NaN;
% grp(11).ID       = {'JPctd_19_1'};
% grp(11).is_2hz   = true;
% grp(11).rsk_file = fans_rsk_JP;
% 
% grp(12).name     = 'JPctd_20';
% grp(12).profiles = 2:4;
% grp(12).lat      = [NaN, NaN, NaN];
% grp(12).lon      = [NaN, NaN, NaN];
% grp(12).ID       = {'JPctd_20_1','JPctd_20_2','JPctd_20_3'};
% grp(12).is_2hz   = true;
% grp(12).rsk_file = fans_rsk_JP;
% 
% % ICM test station (16 Hz)
% grp(13).name     = 'ICMctd_20';
% grp(13).profiles = 13:14;  % v3 ref: 13:14 (intercalated between 2C and 3D) — v4 had 43:44 (WRONG)
% grp(13).lat      = [NaN, NaN];
% grp(13).lon      = [NaN, NaN];
% grp(13).ID       = {'ICM20_1','ICM20_2'};
% grp(13).is_2hz   = false;
% grp(13).rsk_file = fans_rsk_ICM;

% Fjord names per group
grp(1).fjord  = 'Sisimiut Fjord';          % F1A
grp(2).fjord  = 'Sisimiut Fjord';          % F1B
grp(3).fjord  = 'Sisimiut Fjord';          % F2C
grp(4).fjord  = 'Ata Fjord';               % F3D
grp(5).fjord  = 'Narsarsuaaraq Fjord';     % F4E
grp(6).fjord  = 'Narsarsuaaraq Fjord';     % F5F
grp(7).fjord  = 'Qaqortoq Fjord';          % F6G
grp(8).fjord  = 'Aappilattoq Fjord';       % F7H
grp(9).fjord  = 'Aappilattoq Fjord';       % F7I
grp(10).fjord = 'Prince Christian Glacier'; % F8J

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
last_rsk_file = '';
rsk = [];

for ig = 1:length(grp)
    g = grp(ig);
    fprintf('\n--- Group %s ---\n', g.name)

    % Open the .rsk file only if it changes
    % (avoids reopening the large ICM file for every group)
    if ~strcmp(g.rsk_file, last_rsk_file)
        if exist('rsk', 'var')
            clear rsk   % close previous SQLite connection
        end
        fprintf('  Opening: %s\n', g.rsk_file)
        rsk = RSKopen(g.rsk_file);
        last_rsk_file = g.rsk_file;
    end

    % Read only the profiles for this group (downcast)
    rsk_g = RSKreadprofiles(rsk, 'profile', g.profiles, 'direction', 'down');
    np_g = length(rsk_g.data);
    fprintf('  %d profiles read (RSK indices: %d:%d)\n', np_g, g.profiles(1), g.profiles(end))

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
end

if exist('rsk', 'var')
    clear rsk   % close final SQLite connection
end

fprintf('\n  Total: %d profiles extracted\n', ip)
if ip ~= n_total
    warning('Expected %d profiles, got %d', n_total, ip)
end

%% ════════════════════════════════════════════════════════════════════════
%  SAVE
%  ════════════════════════════════════════════════════════════════════════

outfile = [out_path 'RAW_CTD_FANS_oceanCasts.mat'];
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
disp('║  FANS RAW — DONE                                                   ║')
disp('╚════════════════════════════════════════════════════════════════════╝')
summary_parts = arrayfun(@(g) sprintf('%s(%d)', g.name, length(g.ID)), grp, 'UniformOutput', false);
fprintf('  %d profiles | %s\n', n_total, strjoin(summary_parts, ' + '))
fprintf('  n_samp min/max: %d / %d\n', min(n_samp), max(n_samp))
fprintf('  Frequencies: 16 Hz (ICM) | 2 Hz (JP)\n')
