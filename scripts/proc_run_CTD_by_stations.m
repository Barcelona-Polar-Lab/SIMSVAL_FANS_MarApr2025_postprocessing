%% proc_run_CTD_by_stations.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nina HOAREAU, February 2026 (BPL, ICM-CSIC)
% ---
% Entry point for per-station L2 CTD processing with integrated soak
% validation. Run this script and answer the interactive prompts to
% configure the campaign, station, soak method, and binning option.
%
% Calls: process_RBR_CTD_v5()
% Output: PROC_CTD_<CAMPAIGN>_<STATION>.mat in ../outputs/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% Toolboxes — Mac paths (active)
% addpath('/Users/ninah/Desktop/WORK/MATLAB_toolbox/rbr-rsktools/');
% addpath('/Users/ninah/Desktop/WORK/MATLAB_toolbox/mksqlite/');
% addpath('/Users/ninah/Desktop/WORK/MATLAB_toolbox/teos10/');
% addpath('/Users/ninah/Desktop/WORK/MATLAB_toolbox/teos10/library/');
% addpath('/Users/ninah/Desktop/WORK/MATLAB_toolbox/teos10/thermodynamics_from_t/');
% addpath('/Users/ninah/Desktop/WORK/MATLAB_toolbox/cmocean/');

% Linux paths (comment out Mac block above and uncomment below)
addpath('/home/nina/Escritorio/WORK/Programs/MATLAB_toolbox/rbr-rsktools/');
addpath('/home/nina/Escritorio/WORK/Programs/MATLAB_toolbox/mksqlite/');
addpath('/home/nina/Escritorio/WORK/Programs/MATLAB_toolbox/teos10/');
addpath('/home/nina/Escritorio/WORK/Programs/MATLAB_toolbox/teos10/library/');
addpath('/home/nina/Escritorio/WORK/Programs/MATLAB_toolbox/teos10/thermodynamics_from_t/');
addpath('/home/nina/Escritorio/WORK/Programs/MATLAB_toolbox/cmocean/');

%% ════════════════════════════════════════════════════════════════════════
%  STATION DEFINITIONS  (do not edit unless coordinates change)
%  ════════════════════════════════════════════════════════════════════════

% --- SIMSVAL ---
simsval_stations.A.lat = [60.699583, 60.695667, 60.688550, 60.693117, ...
                          60.695383, 60.694983, 60.692967, 60.673933];
simsval_stations.A.lon = [-46.082883, -46.075933, -46.075333, -46.094367, ...
                          -46.106333, -46.160033, -46.165383, -46.156833];
simsval_stations.A.ID = {'A1','A2','A3','A4','A5','A6','A7','A8'};
simsval_stations.A.rsk = '../RAW/raw_SIMSVAL/20250324_A.rsk';
simsval_stations.A.atm = '../ancillary_data/Atmospheric_Data_SIMSVAL_hourlyMean.mat';
simsval_stations.A.comment = 'Qaqortoq Fjord';
simsval_stations.A.fjord  = 'Qaqortoq Fjord';
simsval_stations.A.is_2hz = false;
simsval_stations.A.sss = [32.49428, 32.52720, 32.66762, 32.51477, 32.51953, NaN, 32.51924, 32.66707];

simsval_stations.B.lat = [60.779000, 60.786750, 60.796450];
simsval_stations.B.lon = [-45.629200, -45.634583, -45.639683];
simsval_stations.B.ID = {'B1','B2','B3'};
simsval_stations.B.rsk = '../RAW/raw_SIMSVAL/20250324_B.rsk';
simsval_stations.B.atm = '../ancillary_data/Atmospheric_Data_SIMSVAL_hourlyMean.mat';
simsval_stations.B.comment = 'Eqaluit Fjord';
simsval_stations.B.fjord  = 'Eqaluit Fjord';
simsval_stations.B.is_2hz = false;
simsval_stations.B.sss = [32.51156, 32.52050, 32.51875];

simsval_stations.B1_val.lat = 60.779000;
simsval_stations.B1_val.lon = -45.629200;
simsval_stations.B1_val.ID = {'B1_val'};
simsval_stations.B1_val.rsk = '../RAW/raw_SIMSVAL/20250324_B1_val.rsk';
simsval_stations.B1_val.atm = '../ancillary_data/Atmospheric_Data_SIMSVAL_hourlyMean.mat';
simsval_stations.B1_val.comment = 'Eqaluit Fjord';
simsval_stations.B1_val.fjord  = 'Eqaluit Fjord';
simsval_stations.B1_val.is_2hz = true;
simsval_stations.B1_val.sss = NaN;

simsval_stations.C.lat = [60.294683, 60.298683, 60.302633, 60.346950, 60.346267, ...
                          60.345617, 60.265150, 60.241700, 60.229583, 60.228817, 60.224400];
simsval_stations.C.lon = [-44.004317, -44.006933, -44.008333, -44.064733, -44.078833, ...
                          -44.093250, -44.094050, -44.131533, -44.163150, -44.157800, -44.145083];
simsval_stations.C.ID = {'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11'};
simsval_stations.C.rsk = '../RAW/raw_SIMSVAL/20250326_C.rsk';
simsval_stations.C.atm = '../ancillary_data/Atmospheric_Data_SIMSVAL_hourlyMean.mat';
simsval_stations.C.comment = 'Aappilattoq Fjord';
simsval_stations.C.fjord  = 'Aappilattoq Fjord';
simsval_stations.C.is_2hz = true;
simsval_stations.C.sss = [32.81570, 32.81278, 32.79778, 32.60832, 32.71481, 32.58807, 32.81611, NaN, 32.78207, NaN, 32.80570];

simsval_stations.I.lat = [65.985000, 65.986517, 65.992333, 65.984300, ...
                          65.948283, 65.906017, 65.912317, 65.920100, 65.943533];
simsval_stations.I.lon = [-52.523333, -52.530700, -52.537200, -52.590533, ...
                          -52.644367, -52.816067, -52.822633, -52.825783, -52.722583];
simsval_stations.I.ID = {'I1','I2','I3','I4','I5','I6','I7','I8','I9'};
simsval_stations.I.rsk = '../RAW/raw_SIMSVAL/237957_20250404_I.rsk';
simsval_stations.I.atm = '../ancillary_data/Atmospheric_Data_SIMSVAL_hourlyMean.mat';
simsval_stations.I.comment = 'Eternity Fjord';
simsval_stations.I.fjord  = 'Eternity Fjord';
simsval_stations.I.is_2hz = false;
simsval_stations.I.sss = [32.62136, 32.61552, 32.60579, NaN, NaN, 32.73082, 32.65460, 32.64818, NaN];

% Sea ice stations (16 Hz)
simsval_stations.D.lat = [60.101, 60.099];
simsval_stations.D.lon = [-43.423, -43.425];
simsval_stations.D.ID = {'D1','D2'};
simsval_stations.D.rsk = '../RAW/raw_SIMSVAL/237957_20250326_D.rsk';
simsval_stations.D.atm = '../ancillary_data/Atmospheric_Data_SIMSVAL_hourlyMean.mat';
simsval_stations.D.comment = 'Sea ice station D';
simsval_stations.D.fjord  = '';
simsval_stations.D.is_2hz = false;

simsval_stations.G.lat = [69.500, 69.502];
simsval_stations.G.lon = [-54.100, -54.105];
simsval_stations.G.ID = {'G1','G2'};
simsval_stations.G.rsk = '../RAW/raw_SIMSVAL/237957_20250401_G.rsk';
simsval_stations.G.atm = '../ancillary_data/Atmospheric_Data_SIMSVAL_hourlyMean.mat';
simsval_stations.G.comment = 'Sea ice station G';
simsval_stations.G.fjord  = '';
simsval_stations.G.is_2hz = false;

simsval_stations.H.lat = [68.900, 68.902];
simsval_stations.H.lon = [-59.000, -59.005];
simsval_stations.H.ID = {'H1','H2'};
simsval_stations.H.rsk = '../RAW/raw_SIMSVAL/237957_20250402_H.rsk';
simsval_stations.H.atm = '../ancillary_data/Atmospheric_Data_SIMSVAL_hourlyMean.mat';
simsval_stations.H.comment = 'Sea ice station H';
simsval_stations.H.fjord  = '';
simsval_stations.H.is_2hz = false;

% --- FANS ---
fans_rsk = '../RAW/raw_FANS/237957_20250429_1235.rsk';
fans_atm = '../ancillary_data/Atmospheric_Data_FANS_hourlyMean.mat';

fans_stations.F1A.lat = [66.9417583, 66.9428383, 66.9447300, 66.9447000, 66.9437967];
fans_stations.F1A.lon = [-53.7110900, -53.7285250, -53.7442583, -53.7461067, -53.7461067];
fans_stations.F1A.ID = {'1A_1','1A_2','1A_3','1A_4','1A_5'};
fans_stations.F1A.profiles = 2:6;   % profile 1 = failed cast A1 (v3 ref)
fans_stations.F1A.comment = 'Sisimiut Fjord';
fans_stations.F1A.fjord  = 'Sisimiut Fjord';
fans_stations.F1A.is_2hz = false;
fans_stations.F1A.sss = [32.87711, 32.91798, 32.82915, 32.83352, 32.90435];

fans_stations.F1B.lat = [66.9419883, 66.9417800, 66.9574667];
fans_stations.F1B.lon = [-53.7690100, -53.7628567, -53.7525000];
fans_stations.F1B.ID = {'1B_1','1B_2','1B_3'};
fans_stations.F1B.profiles = 7:9;
fans_stations.F1B.comment = 'Sisimiut Fjord';
fans_stations.F1B.fjord  = 'Sisimiut Fjord';
fans_stations.F1B.is_2hz = false;
fans_stations.F1B.sss = [32.81862, 32.84053, 32.82049];

fans_stations.F2C.lat = [66.8891983, 66.8931433, 66.8789650];
fans_stations.F2C.lon = [-53.8549033, -53.6419133, -53.6877367];
fans_stations.F2C.ID = {'2C_1','2C_2','2C_3'};
fans_stations.F2C.profiles = 10:12;
fans_stations.F2C.comment = 'Sisimiut Fjord';
fans_stations.F2C.fjord  = 'Sisimiut Fjord';
fans_stations.F2C.is_2hz = false;
fans_stations.F2C.sss = [32.81142, NaN, 32.82837];

fans_stations.F3D.lat = [69.4341750, 69.4532033, 69.4938617, 69.4328983];
fans_stations.F3D.lon = [-50.8942933, -51.0013517, -51.0517317, -51.0828333];
fans_stations.F3D.ID = {'3D_1','3D_3','3D_4','3D_5'};
fans_stations.F3D.profiles = [15, 17, 18, 19];  % profiles 13-14 = ICMctd_20 intercalated (v3 ref); profile 16 = 3D_2 excluded (failed cast, unknown position)
fans_stations.F3D.comment = 'Ata Fjord';
fans_stations.F3D.fjord  = 'Ata Fjord';
fans_stations.F3D.is_2hz = false;
fans_stations.F3D.sss = [32.54020, NaN, NaN, NaN];

fans_stations.F4E.lat = [61.1372833, 61.1269733, 61.1270667, 61.1283317, 61.1291317];
fans_stations.F4E.lon = [-45.4901250, -45.4807433, -45.4651233, -45.4426867, -45.4210700];
fans_stations.F4E.ID = {'4E_1','4E_2','4E_3','4E_4','4E_5'};
fans_stations.F4E.profiles = 21:25;  % profile 20 = transit cast between 3D and 4E (v3 ref)
fans_stations.F4E.comment = 'Narsarsuaaraq Fjord';
fans_stations.F4E.fjord  = 'Narsarsuaaraq Fjord';
fans_stations.F4E.is_2hz = false;
fans_stations.F4E.sss = [NaN, NaN, NaN, NaN, NaN];

fans_stations.F5F.lat = [61.0162367, 61.0229967, 61.0316933, 61.0431567, 61.0040883];
fans_stations.F5F.lon = [-46.1302633, -46.1286933, -46.1269300, -46.1303733, -46.0789850];
fans_stations.F5F.ID = {'5F_1','5F_2','5F_3','5F_4','5F_5'};
fans_stations.F5F.profiles = 26:30;
fans_stations.F5F.comment = 'Narsarsuaaraq Fjord';
fans_stations.F5F.fjord  = 'Narsarsuaaraq Fjord';
fans_stations.F5F.is_2hz = false;
fans_stations.F5F.sss = [25.01646, NaN, NaN, NaN, NaN];

fans_stations.F6G.lat = [60.7106683, 60.7056267, 60.7000133, 60.6918383];
fans_stations.F6G.lon = [-46.0123300, -46.0086967, -46.0044550, -46.0022920];
fans_stations.F6G.ID = {'6G_1','6G_2','6G_3','6G_4'};
fans_stations.F6G.profiles = 31:34;
fans_stations.F6G.comment = 'Qaqortoq Fjord';
fans_stations.F6G.fjord  = 'Qaqortoq Fjord';
fans_stations.F6G.is_2hz = false;
fans_stations.F6G.sss = [NaN, NaN, NaN, NaN];

fans_stations.F7H.lat = [60.2873383, 60.2913917, 60.2945233, 60.2967983];
fans_stations.F7H.lon = [-44.2383067, -44.2293433, -44.2200617, -44.2125383];
fans_stations.F7H.ID = {'7H_1','7H_2','7H_3','7H_4'};
fans_stations.F7H.profiles = 35:38;
fans_stations.F7H.comment = 'Aappilattoq Fjord';
fans_stations.F7H.fjord  = 'Aappilattoq Fjord';
fans_stations.F7H.is_2hz = false;
fans_stations.F7H.sss = [NaN, NaN, NaN, NaN];

fans_stations.F7I.lat = [60.2638117, 60.2625383, 60.2603300];
fans_stations.F7I.lon = [-44.1738883, -44.1916250, -44.2118333];
fans_stations.F7I.ID = {'7I_1','7I_2','7I_3'};
fans_stations.F7I.profiles = 39:41;
fans_stations.F7I.comment = 'Aappilattoq Fjord';
fans_stations.F7I.fjord  = 'Aappilattoq Fjord';
fans_stations.F7I.is_2hz = false;
fans_stations.F7I.sss = [NaN, NaN, NaN];

fans_stations.F8J.lat = [60.1801050, 60.1791633, 60.1782383, 60.1801383, 60.1609800];
fans_stations.F8J.lon = [-43.6296067, -43.6093467, -43.6203117, -43.6380100, -43.6284800];
fans_stations.F8J.ID = {'8J_1','8J_2','8J_3','8J_4','8J_5'};
fans_stations.F8J.profiles = 42:46;
fans_stations.F8J.comment = 'Prince Christian Glacier';
fans_stations.F8J.fjord  = 'Prince Christian Glacier';
fans_stations.F8J.is_2hz = false;
fans_stations.F8J.sss = [29.53543, NaN, NaN, NaN, NaN];

%% ════════════════════════════════════════════════════════════════════════
%  INTERACTIVE CONFIGURATION
%  ════════════════════════════════════════════════════════════════════════

fprintf('\n');
disp('╔════════════════════════════════════════════════════════════════════╗');
disp('║  CTD PROCESSING v5 — Interactive configuration                    ║');
disp('╚════════════════════════════════════════════════════════════════════╝');
fprintf('\n');

% ── [1/5]  Campaign ──────────────────────────────────────────────────────
disp('  ── [1/5]  Campaign ──────────────────────────────────────────────────');
disp('  Available: SIMSVAL  FANS');
while true
    campaign_type = upper(strtrim(input('  Campaign [SIMSVAL]: ', 's')));
    if isempty(campaign_type), campaign_type = 'SIMSVAL'; end
    if ismember(campaign_type, {'SIMSVAL', 'FANS'}), break; end
    fprintf('  ✗  Enter SIMSVAL or FANS.\n');
end

% ── [2/5]  Station ───────────────────────────────────────────────────────
fprintf('\n');
disp('  ── [2/5]  Station ───────────────────────────────────────────────────');
if strcmp(campaign_type, 'SIMSVAL')
    avail_list = fieldnames(simsval_stations);
    avail_str  = strjoin(avail_list, '  ');
else
    avail_list = regexprep(fieldnames(fans_stations), '^F', '');
    avail_str  = strjoin(avail_list, '  ');
end
fprintf('  Available: %s\n', avail_str);
while true
    station = strtrim(input('  Station: ', 's'));
    if strcmp(campaign_type, 'SIMSVAL')
        if isfield(simsval_stations, station), break; end
    else
        fn = station;
        if ~startsWith(fn,'F') && ~startsWith(fn,'ICM') && ~startsWith(fn,'JP')
            fn = ['F' fn];
        end
        if isfield(fans_stations, fn), break; end
    end
    fprintf('  ✗  Unknown station. Available: %s\n', avail_str);
end

% ── [3/5]  Soak detection ────────────────────────────────────────────────
fprintf('\n');
disp('  ── [3/5]  Soak detection ────────────────────────────────────────────');
disp('    fixed_time  — remove first N seconds        (recommended: FANS)');
disp('    velocity    — automatic via descent velocity (recommended: SIMSVAL)');
while true
    soak_method = lower(strtrim(input('  Method [fixed_time]: ', 's')));
    if isempty(soak_method), soak_method = 'fixed_time'; end
    if ismember(soak_method, {'fixed_time', 'velocity'}), break; end
    fprintf('  ✗  Enter fixed_time or velocity.\n');
end

if strcmp(soak_method, 'fixed_time')
    val = strtrim(input('  fixed_time_s [20 s]: ', 's'));
    if isempty(val), soak_fixed_time_s = 20; else, soak_fixed_time_s = str2double(val); end
else
    val = strtrim(input('  velocity_threshold [0.15 m/s]: ', 's'));
    if isempty(val), soak_vel_threshold = 0.15; else, soak_vel_threshold = str2double(val); end
    fprintf('  (window_size and min_consecutive are auto-adapted to sampling frequency)\n');
end

val = strtrim(input('  min_soak_depth [0.0 dbar]: ', 's'));
if isempty(val), min_soak_depth = 0.0; else, min_soak_depth = str2double(val); end

% ── [4/5]  Binning ───────────────────────────────────────────────────────
fprintf('\n');
disp('  ── [4/5]  Binning ───────────────────────────────────────────────────');
val = lower(strtrim(input('  Bin to 0.25 dbar grid? [y/N]: ', 's')));
do_binning = strcmp(val, 'y');
if do_binning
    val = strtrim(input('  Bin size [0.25 dbar]: ', 's'));
    if isempty(val), binning_size = 0.25; else, binning_size = str2double(val); end
    val = strtrim(input('  Pressure range [0 150] dbar  (e.g. "0 200"): ', 's'));
    if isempty(val)
        binning_boundary = [0 150];
    else
        binning_boundary = str2num(val); %#ok<ST2NM>
    end
end

% ── [5/5]  Options ───────────────────────────────────────────────────────
fprintf('\n');
disp('  ── [5/5]  Options ───────────────────────────────────────────────────');
val = lower(strtrim(input('  Save soak validation figures? [Y/n]: ', 's')));
save_figures = ~strcmp(val, 'n');

val = lower(strtrim(input('  Interactive mode (pause on suspect profiles)? [y/N]: ', 's')));
interactive_mode = strcmp(val, 'y');

soak_thresholds.depth_max = 5;  % suspect threshold — edit in script if needed

% ── Summary ──────────────────────────────────────────────────────────────
fprintf('\n');
disp('  ── Summary ──────────────────────────────────────────────────────────');
fprintf('  Campaign  : %s\n', campaign_type);
fprintf('  Station   : %s\n', station);
if strcmp(soak_method, 'fixed_time')
    fprintf('  Soak      : fixed_time  (%.0f s,  min_depth %.1f dbar)\n', soak_fixed_time_s, min_soak_depth);
else
    fprintf('  Soak      : velocity    (threshold %.2f m/s,  min_depth %.1f dbar)\n', soak_vel_threshold, min_soak_depth);
end
if do_binning
    fprintf('  Binning   : %.2f dbar  [%.0f–%.0f dbar]\n', binning_size, binning_boundary(1), binning_boundary(2));
else
    fprintf('  Binning   : native resolution\n');
end
fprintf('  Figures   : %s\n', iif(save_figures, 'yes', 'no'));
fprintf('  Mode      : %s\n', iif(interactive_mode, 'interactive', 'automatic'));
fprintf('\n');
input('  Press Enter to start, or Ctrl+C to abort...', 's');
fprintf('\n');

% ── Build opts struct ─────────────────────────────────────────────────────
opts.cast            = 'down';
opts.do_binning      = do_binning;
if do_binning
    opts.binning     = binning_size;
    opts.boundary    = binning_boundary;
end
opts.interactive     = interactive_mode;
opts.save_figures    = save_figures;
opts.soak_thresholds = soak_thresholds;
opts.min_soak_depth  = min_soak_depth;
opts.soak_method     = soak_method;
if strcmp(soak_method, 'fixed_time')
    opts.soak_fixed_time_s  = soak_fixed_time_s;
else
    opts.soak_vel_threshold = soak_vel_threshold;
end

% --- Per-profile soak overrides (optional) ---
% Override soak parameters for specific profiles. Leave empty for none.
% Supported fields per profile: .method, .fixed_time_s, .vel_threshold
% Example — profile 3 needs 40s instead of default, profile 5 uses velocity:
%   opts.soak_overrides(3).fixed_time_s = 40;
%   opts.soak_overrides(5).method       = 'velocity';
%   opts.soak_overrides(5).vel_threshold = 0.10;
opts.soak_overrides = [];

%% ════════════════════════════════════════════════════════════════════════
%  PARAMETER PREPARATION
%  ════════════════════════════════════════════════════════════════════════

switch upper(campaign_type)
    case 'SIMSVAL'
        st = simsval_stations.(station);
        station_info.lat      = st.lat;
        station_info.lon      = st.lon;
        station_info.ID       = st.ID;
        station_info.campaign = ['SIMSVAL_' station];
        station_info.vessel   = 'Le Commandant Charcot';
        station_info.comment  = st.comment;
        rsk_file          = st.rsk;
        opts.atm_file     = st.atm;
        opts.output_path  = '../outputs/';
        opts.fig_path     = '../outputs/validation_soak_figures/';
        opts.is_2hz       = st.is_2hz;
        if isfield(st, 'sss'),   station_info.sss   = st.sss;   end
        if isfield(st, 'fjord'), station_info.fjord = st.fjord; end

    case 'FANS'
        fn = station;
        if ~startsWith(fn,'F') && ~startsWith(fn,'ICM') && ~startsWith(fn,'JP')
            fn = ['F' fn];
        end
        st = fans_stations.(fn);
        station_info.lat      = st.lat;
        station_info.lon      = st.lon;
        station_info.ID       = st.ID;
        station_info.campaign = ['FANS_' station];
        station_info.vessel   = 'Le Commandant Charcot';
        station_info.comment  = st.comment;
        rsk_file = fans_rsk;
        if isfield(st, 'rsk'), rsk_file = st.rsk; end
        opts.profiles    = st.profiles;
        opts.atm_file    = fans_atm;
        opts.output_path = '../outputs/';
        opts.fig_path    = '../outputs/validation_soak_figures/';
        opts.is_2hz      = st.is_2hz;
        if isfield(st, 'sss'),   station_info.sss   = st.sss;   end
        if isfield(st, 'fjord'), station_info.fjord = st.fjord; end
end

%% ════════════════════════════════════════════════════════════════════════
%  PER-PROFILE SOAK OVERRIDES (optional interactive step)
%  ════════════════════════════════════════════════════════════════════════

fprintf('\n');
disp('─────────────────────────────────────────────────────────────────');
disp('  [Optional] Per-profile soak overrides');
fprintf('  Default: %s', opts.soak_method);
if strcmp(opts.soak_method, 'fixed_time')
    fprintf(' (%.0f s)', opts.soak_fixed_time_s);
else
    fprintf(' (threshold %.2f m/s)', opts.soak_vel_threshold);
end
fprintf('\n');
disp('  Available profiles:');
for k = 1:length(station_info.ID)
    fprintf('    [%d]  %s\n', k, station_info.ID{k});
end
disp('─────────────────────────────────────────────────────────────────');

ans_ov = strtrim(input('  Override any profile? [no]: ', 's'));
if ismember(lower(ans_ov), {'y', 'yes', 'oui', 'o'})
    while true
        prof_str = strtrim(input('  Profile number (Enter to finish): ', 's'));
        if isempty(prof_str), break; end
        prof_idx = str2double(prof_str);
        if isnan(prof_idx) || prof_idx < 1 || prof_idx > length(station_info.ID)
            fprintf('  ✗  Invalid profile number. Enter a number between 1 and %d.\n', length(station_info.ID));
            continue;
        end
        prof_idx = round(prof_idx);
        fprintf('  Profile [%d] %s — current: %s', prof_idx, station_info.ID{prof_idx}, opts.soak_method);
        if strcmp(opts.soak_method, 'fixed_time')
            fprintf(' (%.0f s)', opts.soak_fixed_time_s);
        end
        fprintf('\n');

        % Method
        while true
            new_method = lower(strtrim(input(sprintf('    Method [%s]: ', opts.soak_method), 's')));
            if isempty(new_method), new_method = opts.soak_method; end
            if ismember(new_method, {'fixed_time', 'velocity'}), break; end
            fprintf('    ✗  Enter fixed_time or velocity.\n');
        end
        opts.soak_overrides(prof_idx).method = new_method;

        % Method-specific parameter
        if strcmp(new_method, 'fixed_time')
            default_ft = opts.soak_fixed_time_s;
            val = strtrim(input(sprintf('    fixed_time_s [%.0f s]: ', default_ft), 's'));
            if isempty(val)
                opts.soak_overrides(prof_idx).fixed_time_s = default_ft;
            else
                opts.soak_overrides(prof_idx).fixed_time_s = str2double(val);
            end
            fprintf('  ✓  Profile [%d] %s → fixed_time_s = %.0f s\n', ...
                prof_idx, station_info.ID{prof_idx}, opts.soak_overrides(prof_idx).fixed_time_s);
        else
            default_vt = opts.soak_vel_threshold;
            val = strtrim(input(sprintf('    vel_threshold [%.2f m/s]: ', default_vt), 's'));
            if isempty(val)
                opts.soak_overrides(prof_idx).vel_threshold = default_vt;
            else
                opts.soak_overrides(prof_idx).vel_threshold = str2double(val);
            end
            fprintf('  ✓  Profile [%d] %s → velocity threshold = %.2f m/s\n', ...
                prof_idx, station_info.ID{prof_idx}, opts.soak_overrides(prof_idx).vel_threshold);
        end
    end
    n_ov = sum(arrayfun(@(x) ~isempty(x), opts.soak_overrides));
    fprintf('\n  %d override(s) defined.\n', n_ov);
else
    disp('  No overrides — default soak applied to all profiles.');
end

%% ════════════════════════════════════════════════════════════════════════
%  RUN PROCESSING
%  ════════════════════════════════════════════════════════════════════════

if ~exist(rsk_file, 'file')
    error('RSK file not found: %s', rsk_file);
end

output = process_RBR_CTD(rsk_file, station_info, opts);

%% ════════════════════════════════════════════════════════════════════════
%  DONE
%  ════════════════════════════════════════════════════════════════════════

fprintf('\n');
disp('╔════════════════════════════════════════════════════════════════════╗');
disp('║  PROCESSING COMPLETE                                              ║');
disp('╚════════════════════════════════════════════════════════════════════╝');
fprintf('\n');
fprintf('  Output file     : %s\n', output.matfile);
fprintf('  Profiles done   : %d / %d\n', output.np_processed, output.np_total);
if output.np_skipped > 0
    fprintf('  Profiles skipped: %d\n', output.np_skipped);
end
fprintf('\n');
disp('  To process another station, run this script again.');
fprintf('\n');

%% ────────────────────────────────────────────────────────────────────────
% Helper: inline if  (avoids anonymous function overhead for display only)
function out = iif(cond, a, b)
    if cond, out = a; else, out = b; end
end
