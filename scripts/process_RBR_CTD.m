function [output] = process_RBR_CTD_v5(rsk_file, station_info, opts)
%% PROCESS_RBR_CTD — RBR CTD processing with integrated soak validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nina HOAREAU, February 2026 (BPL, ICM-CSIC)
% ---
% Processes a single RBR Concerto .rsk file through the full L2 pipeline:
%   RAW read → sea pressure → soak removal → corrections → QC flags → export
%
% USAGE:
%   output = process_RBR_CTD_v5(rsk_file, station_info, opts)
%
% INPUTS:
%   rsk_file     - Path to the .rsk file
%   station_info - Structure with:
%       .lat      - Latitude vector (1 × np)
%       .lon      - Longitude vector (1 × np)
%       .ID       - Cell array of station IDs
%       .campaign - Campaign name (e.g. 'SIMSVAL_A')
%       .vessel   - Vessel name
%       .comment  - Comment string
%   opts - Options structure:
%       .cast           - 'down' or 'up'             (default: 'down')
%       .patm           - Atmospheric pressure (dbar) (default: 10.1325)
%       .atm_file       - Atmospheric data file       (optional)
%       .binning        - Bin size (dbar)             (default: 0.25)
%       .boundary       - Pressure bounds [min max]  (default: [0 150])
%       .profiles       - Profiles to process        (default: all)
%       .do_binning     - Apply RSKbinaverage         (default: true)
%       .interactive    - Interactive mode            (default: false)
%       .save_figures   - Save validation figures     (default: true)
%       .fig_path       - Figures output path        (default: './validation_figures/')
%       .output_path    - .mat output path           (default: './')
%       .soak_thresholds - Suspect soak thresholds:
%           .depth_max   - Max soak depth (dbar)     (default: 10)
%
% OUTPUTS:
%   output - Structure with:
%       .matfile        - Generated .mat file path
%       .np_total       - Total number of profiles
%       .np_processed   - Number of processed profiles
%       .np_skipped     - Number of skipped profiles
%       .validation     - Per-profile validation results
%       .qc_summary     - QC summary
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0. Initialisation and default parameters
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%  PROCESS_RBR_CTD — CTD processing with integrated soak validation    %')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')

% Default parameters
if ~isfield(opts, 'cast'),        opts.cast        = 'down';   end
if ~isfield(opts, 'patm'),        opts.patm        = 10.1325;  end
if ~isfield(opts, 'binning'),     opts.binning     = 0.25;     end
if ~isfield(opts, 'boundary'),    opts.boundary    = [0 150];  end
if ~isfield(opts, 'profiles'),    opts.profiles    = [];       end
if ~isfield(opts, 'do_binning'),  opts.do_binning  = true;     end
if ~isfield(opts, 'interactive'), opts.interactive = false;    end
if ~isfield(opts, 'save_figures'),opts.save_figures = true;    end
if ~isfield(opts, 'fig_path'),    opts.fig_path    = './validation_figures/'; end
if ~isfield(opts, 'output_path'), opts.output_path = './';     end

% Default soak thresholds
if ~isfield(opts, 'soak_thresholds')
    opts.soak_thresholds.depth_max = 10;
end

% Velocity thresholds (propagated to all figures)
if ~isfield(opts, 'velocity_threshold_soak'), opts.velocity_threshold_soak = 0.15; end  % soak detection (m/s)
if ~isfield(opts, 'velocity_threshold_loop'), opts.velocity_threshold_loop = 0.1;  end  % RSKremoveloops (m/s)

% Create figures directory if needed
if opts.save_figures && ~exist(opts.fig_path, 'dir')
    mkdir(opts.fig_path);
end

%% 1. Read RSK file
disp('════════════════════════════════════════════════════════════════════')
disp('  STEP 1: READ RSK FILE')
disp('════════════════════════════════════════════════════════════════════')

rsk = RSKopen(rsk_file);
RSKprintchannels(rsk);

% Detect temperature channel
channel_names = {rsk.channels.longName};
is_temp1 = any(strcmpi(channel_names, 'Temperature1')) && ~any(strcmpi(channel_names, 'Temperature'));
if any(strcmpi(channel_names, 'Temperature'))
    temp_channel = 'Temperature';
elseif is_temp1
    temp_channel = 'Temperature';
    fprintf('  Temperature1 channel detected — will be renamed to Temperature\n');
else
    error('Temperature channel not found');
end
fprintf('  Temperature channel: %s\n', temp_channel);

% Detect sampling frequency
% NOTE: rsk.continuous.samplingPeriod is in MILLISECONDS
sampling_period_ms = rsk.continuous.samplingPeriod;
sampling_period_s  = sampling_period_ms / 1000;
sampling_freq      = 1 / sampling_period_s;
fprintf('  Frequency: %.1f Hz\n', sampling_freq);

% High-frequency flag: < 8 Hz = 2 Hz mode, >= 8 Hz = 16 Hz mode
is_highfreq = sampling_freq >= 8;
if is_highfreq
    fprintf('  Mode: 16 Hz\n');
else
    fprintf('  Mode: 2 Hz\n');
end

% Number of profiles
np_total = length(rsk.profiles.downcast.tend);
fprintf('  Profiles detected: %d\n', np_total);

% Profiles to process
if isempty(opts.profiles)
    profiles_to_read = 1:np_total;
else
    profiles_to_read = opts.profiles;
end
np = length(profiles_to_read);

% Read profiles
rsk = RSKreadprofiles(rsk, 'profile', profiles_to_read, 'direction', opts.cast);

% Rename Temperature1 → Temperature after readprofiles (2 Hz PONANT)
% RSKreadprofiles re-reads the SQLite database and resets channel names
if is_temp1
    ch = {rsk.channels.longName};
    idx_t1 = find(strcmpi(ch, 'Temperature1'));
    if ~isempty(idx_t1)
        rsk.channels(idx_t1).longName = 'Temperature';
        fprintf('  Temperature1 renamed → Temperature\n');
    end
end

%% 2. Station metadata
disp(' ')
disp('════════════════════════════════════════════════════════════════════')
disp('  STEP 2: STATION METADATA')
disp('════════════════════════════════════════════════════════════════════')

% RSKaddstationdata may fail for single-profile files
if np > 1
    rsk = RSKaddstationdata(rsk, 'profile', profiles_to_read, ...
        'latitude', station_info.lat, 'longitude', station_info.lon, ...
        'station', station_info.ID, 'vessel', station_info.vessel, ...
        'comment', station_info.comment);
else
    disp('  -> Single profile: RSKaddstationdata skipped (manual lat/lon storage)');
end

fprintf('  Campaign: %s\n', station_info.campaign);
fprintf('  Stations: %s\n', strjoin(station_info.ID, ', '));

%% 3. Atmospheric pressure correction → sea pressure
disp(' ')
disp('════════════════════════════════════════════════════════════════════')
disp('  STEP 3: SEA PRESSURE DERIVATION')
disp('════════════════════════════════════════════════════════════════════')

% Load atmospheric data if available (FerryBox hourly means)
if isfield(opts, 'atm_file') && exist(opts.atm_file, 'file')
    atm_data = load(opts.atm_file);
    mtimeI = mean(datetime(rsk.profiles.downcast.tstart, 'ConvertFrom', 'datenum'));
    [~, idx] = min(abs(atm_data.hourlyMean.DateTime - mtimeI));
    hPa1 = atm_data.hourlyMean.AtmPress1_hPa(idx);
    hPa2 = atm_data.hourlyMean.AtmPress2_hPa(idx);
    opts.patm = ((hPa1 + hPa2)/2)/100;
    fprintf('  Patm from file: %.4f dbar\n', opts.patm);

    % Store for per-profile atmospheric export
    atm_info.hPa1 = hPa1;
    atm_info.hPa2 = hPa2;
    atm_info.file = opts.atm_file;
    atm_info.data = atm_data;
else
    fprintf('  Default Patm: %.4f dbar\n', opts.patm);
    atm_info = [];
end

rsk = RSKderiveseapressure(rsk, 'patm', opts.patm);

%% 4. Soak detection and removal
disp(' ')
disp('════════════════════════════════════════════════════════════════════')
disp('  STEP 4: AUTOMATIC SOAK DETECTION')
disp('════════════════════════════════════════════════════════════════════')

% Derive depth and velocity (needed for soak detection)
rsk = RSKderivedepth(rsk);
rsk = RSKderivevelocity(rsk);

% Soak parameters adapted to sampling frequency
% Goal: same durations in seconds regardless of frequency
if is_highfreq  % 16 Hz
    soak_window   = 48;  % ~3s velocity smoothing (reduces jerky descent oscillations)
    soak_consec   = 16;  % ~1s of consecutive points above threshold
    stability_window = 64;  % ~4s for std(T) and std(C) (informational only)
else            % 2 Hz
    soak_window   = 6;   % ~3s velocity smoothing
    soak_consec   = 2;   % ~1s
    stability_window = 8;   % ~4s adapted to 2 Hz
end

fprintf('  Frequency: %.0f Hz\n', sampling_freq);
fprintf('  stability_window: %d pts (%.1f s)\n', stability_window, stability_window/sampling_freq);
fprintf('  soak_consec: %d pts (%.1f s)\n', soak_consec, soak_consec/sampling_freq);

% Apply soak filtering — method and parameters from opts (set in proc_run_CTD_by_stations.m)
if ~isfield(opts, 'min_soak_depth'),     opts.min_soak_depth     = 1.0;          end
if ~isfield(opts, 'soak_method'),        opts.soak_method        = 'fixed_time'; end
if ~isfield(opts, 'soak_fixed_time_s'),  opts.soak_fixed_time_s  = 30;           end
if ~isfield(opts, 'soak_vel_threshold'), opts.soak_vel_threshold = 0.15;         end
if ~isfield(opts, 'soak_overrides'),     opts.soak_overrides     = [];           end

fprintf('  Default method  : %s\n', opts.soak_method);
fprintf('  min_soak_depth  : %.1f dbar\n', opts.min_soak_depth);
if ~isempty(opts.soak_overrides)
    fprintf('  Per-profile overrides defined: %d profile(s)\n', length(opts.soak_overrides));
end

% --- Pre-fetch channel indices needed for manual per-profile trim ---
SPcol = getchannelindex(rsk, 'Sea Pressure');
try
    Vcol = getchannelindex(rsk, 'Velocity');
catch
    Vcol = [];
end

% --- Initialise soak_info output struct ---
soak_info = struct('profile', {}, 'duration_s', {}, 'depth_dbar', {}, ...
                   'n_filtered', {}, 'method', {});

fprintf('  %-8s %-14s %10s %10s %10s\n', 'Profile', 'Method', 'Duration(s)', 'SP(dbar)', 'Pts NaN');
disp('  ─────────────────────────────────────────────────────────');

for ip = 1:np

    % ── Determine soak parameters for this profile ──────────────────────
    has_override = ~isempty(opts.soak_overrides) && ...
                   ip <= numel(opts.soak_overrides) && ...
                   ~isempty(opts.soak_overrides(ip));

    if has_override
        ov = opts.soak_overrides(ip);
        if isfield(ov, 'method')
            ip_method = ov.method;
        else
            ip_method = opts.soak_method;
        end
        if isfield(ov, 'fixed_time_s')
            ip_fixed_time_s = ov.fixed_time_s;
        else
            ip_fixed_time_s = opts.soak_fixed_time_s;
        end
        if isfield(ov, 'vel_threshold')
            ip_vel_threshold = ov.vel_threshold;
        else
            ip_vel_threshold = opts.soak_vel_threshold;
        end
    else
        ip_method        = opts.soak_method;
        ip_fixed_time_s  = opts.soak_fixed_time_s;
        ip_vel_threshold = opts.soak_vel_threshold;
    end

    % ── Compute soak end index ───────────────────────────────────────────
    SP   = rsk.data(ip).values(:, SPcol);
    time = rsk.data(ip).tstamp;
    n_pts = length(SP);
    t_rel = (time - time(1)) * 86400;

    if strcmp(ip_method, 'fixed_time')

        idx_fixed = find(t_rel >= ip_fixed_time_s, 1, 'first');
        if isempty(idx_fixed)
            idx_fixed = n_pts;
            warning('process_RBR_CTD: profile %d — fixed_time_s (%.0fs) exceeds profile duration.', ip, ip_fixed_time_s);
        end
        soak_end_idx = idx_fixed;

    else  % velocity

        if isempty(Vcol)
            error('process_RBR_CTD: Velocity channel not found for velocity soak method.');
        end
        V        = rsk.data(ip).values(:, Vcol);
        V_smooth = movmean(V, soak_window, 'omitnan');
        dt       = (time(2) - time(1)) * 86400;
        skip_pts = max(1, floor(5.0 / dt) + 1);
        search_start = min(skip_pts, n_pts - soak_consec);
        soak_end_idx = search_start;
        for i = search_start:(n_pts - soak_consec)
            if all(V_smooth(i:i+soak_consec-1) > ip_vel_threshold)
                soak_end_idx = i;
                break;
            end
        end

    end

    % ── Apply min_soak_depth ─────────────────────────────────────────────
    if SP(soak_end_idx) < opts.min_soak_depth
        idx_min = find(SP >= opts.min_soak_depth, 1, 'first');
        if ~isempty(idx_min) && idx_min > soak_end_idx
            soak_end_idx = idx_min;
        end
    end

    % ── Apply NaN trim ───────────────────────────────────────────────────
    rsk.data(ip).values(1:soak_end_idx, :) = NaN;

    % ── Store soak_info ──────────────────────────────────────────────────
    soak_info(ip).profile    = ip;
    soak_info(ip).duration_s = t_rel(soak_end_idx);
    soak_info(ip).depth_dbar = SP(soak_end_idx);
    soak_info(ip).n_filtered = soak_end_idx;
    soak_info(ip).method     = ip_method;

    flag_str = '';
    if has_override, flag_str = ' ◀ override'; end
    fprintf('  %-8d %-14s %10.1f %10.2f %10d%s\n', ip, ip_method, ...
        soak_info(ip).duration_s, soak_info(ip).depth_dbar, soak_end_idx, flag_str);

end

disp('  ─────────────────────────────────────────────────────────');

% Soak summary
fprintf('\n  Soak summary:\n');
fprintf('    Mean duration : %.1f s\n',   mean([soak_info.duration_s]));
fprintf('    Mean depth    : %.2f dbar\n', mean([soak_info.depth_dbar]));

%% 4bis. Integrated soak validation
disp(' ')
disp('════════════════════════════════════════════════════════════════════')
disp('  STEP 4bis: INTEGRATED SOAK VALIDATION')
disp('════════════════════════════════════════════════════════════════════')
disp(' ')
disp('  RAW vs POST-SOAK comparison for each profile')
disp('  Suspect profile detection thresholds:')
fprintf('    Max soak depth: %.1f dbar\n', opts.soak_thresholds.depth_max);
disp(' ')

% Reload RAW data for comparison (before any processing)
rsk_raw = RSKopen(rsk_file);
rsk_raw = RSKreadprofiles(rsk_raw, 'profile', profiles_to_read, 'direction', opts.cast);
% Apply Temperature1 → Temperature rename AFTER readprofiles
raw_channels = {rsk_raw.channels.longName};
if any(strcmpi(raw_channels, 'Temperature1')) && ~any(strcmpi(raw_channels, 'Temperature'))
    idx_t1 = find(strcmpi(raw_channels, 'Temperature1'));
    rsk_raw.channels(idx_t1).longName = 'Temperature';
end
rsk_raw = RSKderiveseapressure(rsk_raw, 'patm', opts.patm);
rsk_raw = RSKderivedepth(rsk_raw);
rsk_raw = RSKderivevelocity(rsk_raw);

% Get channel indices
[SPcol] = getchannelindex(rsk_raw, 'Sea Pressure');
[Tcol]  = getchannelindex(rsk_raw, temp_channel);
[Ccol]  = getchannelindex(rsk_raw, 'Conductivity');
[Vcol]  = getchannelindex(rsk_raw, 'Velocity');

% Structure to store validation results
validation          = struct();
validation.profiles = struct();
validation.suspects = [];
validation.skipped  = [];

% Print table header
fprintf('  %-6s %-12s %10s %10s %10s %10s\n', ...
        'Profile', 'ID', 'Soak(s)', 'Soak(dbar)', '% filtered', 'Status');
disp(repmat('-', 1, 70));

for ip = 1:np

    prof_idx = profiles_to_read(ip);

    % Profile ID
    if iscell(station_info.ID)
        prof_id = station_info.ID{ip};
    else
        prof_id = char(station_info.ID(ip));
    end

    % Soak statistics
    soak_dur    = soak_info(ip).duration_s;
    soak_dep    = soak_info(ip).depth_dbar;
    n_filtered  = soak_info(ip).n_filtered;

    % Compute % filtered
    n_total_raw  = length(rsk_raw.data(ip).values(:, SPcol));
    pct_filtered = 100 * n_filtered / n_total_raw;

    % Determine if suspect (one criterion: soak depth > depth_max)
    is_suspect      = false;
    suspect_reasons = {};

    if soak_dep > opts.soak_thresholds.depth_max
        is_suspect = true;
        suspect_reasons{end+1} = sprintf('depth=%.1f dbar', soak_dep);
    end

    % Store validation info
    validation.profiles(ip).index               = ip;
    validation.profiles(ip).id                  = prof_id;
    validation.profiles(ip).soak_duration_s     = soak_dur;
    validation.profiles(ip).soak_depth_dbar     = soak_dep;
    validation.profiles(ip).pct_filtered        = pct_filtered;
    validation.profiles(ip).is_suspect          = is_suspect;
    validation.profiles(ip).suspect_reasons     = suspect_reasons;
    validation.profiles(ip).decision            = 'OK';   % default
    validation.profiles(ip).manual_filtered_idx = [];     % manually filtered points

    if is_suspect
        validation.suspects(end+1) = ip;
        status = 'SUSPECT';
        fprintf('  %-6d %-12s %10.1f %10.2f %10.1f   \x1b[33m%s\x1b[0m\n', ...
                ip, prof_id, soak_dur, soak_dep, pct_filtered, status);
    else
        status = 'OK';
        fprintf('  %-6d %-12s %10.1f %10.2f %10.1f   \x1b[32m%s\x1b[0m\n', ...
                ip, prof_id, soak_dur, soak_dep, pct_filtered, status);
    end

end

disp(repmat('-', 1, 70));
fprintf('  Suspect profiles: %d / %d\n', length(validation.suspects), np);
if ~isempty(validation.skipped)
    fprintf('  Excluded profiles (SKIP): %d\n', length(validation.skipped));
end
disp(' ')

% Determine profiles to continue processing
profiles_to_process = setdiff(1:np, validation.skipped);
np_process = length(profiles_to_process);

if np_process == 0
    error('All profiles were excluded. Processing aborted.');
end

fprintf('  Continuing with %d profiles...\n', np_process);

%% 5. Post-processing (on valid profiles)
disp(' ')
disp('════════════════════════════════════════════════════════════════════')
disp('  STEP 5: POST-PROCESSING')
disp('════════════════════════════════════════════════════════════════════')

rsk_processed = rsk;

% 5.1 A2D zero-order hold correction
disp('  5.1 A2D zero-order hold correction...')
rsk_processed = RSKcorrecthold(rsk_processed, 'action', 'interp');

% 5.2 Despiking (T and C, 4σ threshold, direction='down', NaN action)
disp('  5.2 Despiking...')
rsk_processed = RSKdespike(rsk_processed, 'channel', 'conductivity', ...
                'direction', opts.cast, 'threshold', 4, 'windowLength', 15, 'action', 'nan');
rsk_processed = RSKdespike(rsk_processed, 'channel', temp_channel, ...
                'direction', opts.cast, 'threshold', 4, 'windowLength', 15, 'action', 'nan');

% 5.3 CT alignment (16 Hz only)
% Cap at abs(lag) > 2 scans to avoid spurious cross-correlation peaks in weakly stratified water
if is_highfreq
    disp('  5.3 CT alignment...')
    lag_scans  = RSKcalculateCTlag(rsk_processed, 'seapressurerange', [5 30]);
    mlag_scans = round(median(-lag_scans));
    if abs(mlag_scans) > 2
        mlag_scans = sign(mlag_scans) * 2;
        fprintf('  CT lag capped at %+d scans (|lag| > 2)\n', mlag_scans);
    end
    rsk_processed = RSKalignchannel(rsk_processed, 'channel', temp_channel, 'lag', mlag_scans);
else
    disp('  5.3 CT alignment skipped (low frequency)')
    mlag_scans = 0;
end

% 5.4 Smoothing (window = 5 for T and C)
disp('  5.4 Smoothing...')
rsk_processed = RSKsmooth(rsk_processed, 'channel', {temp_channel, 'Conductivity'}, 'windowLength', 5);

% 5.5 Loop removal (velocity threshold = opts.velocity_threshold_loop)
disp('  5.5 Loop removal...')
rsk_processed = RSKderivedepth(rsk_processed);
rsk_processed = RSKderivevelocity(rsk_processed);
rsk_processed = RSKremoveloops(rsk_processed, 'threshold', opts.velocity_threshold_loop);

% 5.6 Derived variables
disp('  5.6 Derived variables...')
rsk_processed = RSKderivesalinity(rsk_processed);
rsk_processed = RSKderivesigma(rsk_processed);

% 5.7 Binning (optional — controlled by do_binning flag)
if opts.do_binning
    disp('  5.7 Binning 0.25 dbar...')
    rsk_processed = RSKbinaverage(rsk_processed, 'binBy', 'Sea Pressure', ...
        'binSize', opts.binning, 'boundary', opts.boundary, 'direction', opts.cast);
else
    disp('  5.7 Binning SKIPPED (native resolution)')
end

disp('  [OK] Post-processing complete')

%% 6. QC tests
disp(' ')
disp('════════════════════════════════════════════════════════════════════')
disp('  STEP 6: QUALITY CONTROL (QC)')
disp('════════════════════════════════════════════════════════════════════')

% Extract data for QC
[spcol_qc] = getchannelindex(rsk_processed, {'Sea Pressure'});
[ccol_qc, tcol_qc, scol_qc, dacol_qc] = getchannelindex(rsk_processed, ...
    {'Conductivity', temp_channel, 'Salinity', 'Density Anomaly'});

if opts.do_binning
    m_qc = round(1/opts.binning * max(opts.boundary));
else
    % Native resolution: find max sample count for NaN-padding
    ns_arr = zeros(1, np);
    for ip = profiles_to_process
        profind_tmp = getdataindex(rsk_processed, 'direction', opts.cast, 'profile', profiles_to_read(ip));
        ns_arr(ip) = size(rsk_processed.data(profind_tmp).values, 1);
    end
    m_qc = max(ns_arr);
end

qc_data.seapres = NaN(m_qc, np);
qc_data.temp    = NaN(m_qc, np);
qc_data.cond    = NaN(m_qc, np);
qc_data.sal     = NaN(m_qc, np);
qc_data.sigma   = NaN(m_qc, np);
qc_data.lat     = zeros(1, np);
qc_data.lon     = zeros(1, np);

for ip = profiles_to_process
    profind_qc = getdataindex(rsk_processed, 'direction', opts.cast, 'profile', profiles_to_read(ip));
    ns = size(rsk_processed.data(profind_qc).values, 1);
    qc_data.seapres(1:ns, ip) = rsk_processed.data(profind_qc).values(:, spcol_qc);
    qc_data.temp(1:ns, ip)    = rsk_processed.data(profind_qc).values(:, tcol_qc);
    qc_data.cond(1:ns, ip)    = rsk_processed.data(profind_qc).values(:, ccol_qc);
    qc_data.sal(1:ns, ip)     = rsk_processed.data(profind_qc).values(:, scol_qc);
    qc_data.sigma(1:ns, ip)   = rsk_processed.data(profind_qc).values(:, dacol_qc);
    if isfield(rsk_processed.data(profind_qc), 'latitude') && ~isempty(rsk_processed.data(profind_qc).latitude)
        qc_data.lat(ip) = rsk_processed.data(profind_qc).latitude;
        qc_data.lon(ip) = rsk_processed.data(profind_qc).longitude;
    else
        qc_data.lat(ip) = station_info.lat(min(ip, length(station_info.lat)));
        qc_data.lon(ip) = station_info.lon(min(ip, length(station_info.lon)));
    end
end

% Apply QC tests
% flatline_n=32 (~2 s at 16 Hz); skip_gradient for native non-binned data (Option C — QARTOD
% gradient thresholds are designed for 0.25 dbar binned data, not 16 Hz native resolution)
if opts.do_binning
    flatline_n = 5;
    skip_grad  = false;
else
    flatline_n = 32;
    skip_grad  = true;
end

[qc_flags, qc_summary] = apply_QC_tests(qc_data, ...
    'region', 'arctic', ...
    'flatline_n', flatline_n, ...
    'density_tol', 0.03, ...
    'skip_gradient', skip_grad, ...
    'verbose', true);

% Integrate validation decisions into QC flags
for ip = 1:np
    decision = validation.profiles(ip).decision;
    switch decision
        case 'BAD'
            % Set all profile flags to 4 (bad)
            qc_flags.temp_qc(:, ip)  = max(qc_flags.temp_qc(:, ip),  4);
            qc_flags.sal_qc(:, ip)   = max(qc_flags.sal_qc(:, ip),   4);
            qc_flags.cond_qc(:, ip)  = max(qc_flags.cond_qc(:, ip),  4);
            qc_flags.sigma_qc(:, ip) = max(qc_flags.sigma_qc(:, ip), 4);
        case 'SUSPECT'
            % Set surface flags (< soak_depth) to 3 (suspect)
            soak_dep     = validation.profiles(ip).soak_depth_dbar;
            surface_idx  = qc_data.seapres(:, ip) < soak_dep;
            qc_flags.temp_qc(surface_idx, ip)  = max(qc_flags.temp_qc(surface_idx, ip),  3);
            qc_flags.sal_qc(surface_idx, ip)   = max(qc_flags.sal_qc(surface_idx, ip),   3);
            qc_flags.cond_qc(surface_idx, ip)  = max(qc_flags.cond_qc(surface_idx, ip),  3);
    end
end

%% 7. Data export
disp(' ')
disp('════════════════════════════════════════════════════════════════════')
disp('  STEP 7: DATA EXPORT')
disp('════════════════════════════════════════════════════════════════════')

m       = m_qc;
ctddate = NaT(1, np);
ctdid   = strings(1, np);
ctdlat  = zeros(1, np);
ctdlon  = zeros(1, np);
ctdcond    = NaN(m, np);
ctdtemp    = NaN(m, np);
ctdsal     = NaN(m, np);
ctdsigma   = NaN(m, np);
ctdvelprof = NaN(m, np);
n_samples  = zeros(1, np);  % number of valid samples per profile

% ctdseapres: 1D if binning (common pressure grid), 2D NaN-padded if native
if opts.do_binning
    ctdseapres = zeros(m, 1);
else
    ctdseapres   = NaN(m, np);
    ctdsampltime = NaN(m, np);  % per-sample timestamps (datenum)
end

% Soak metadata
soak_duration_s = zeros(1, np);
soak_depth_dbar = zeros(1, np);
soak_n_filtered = zeros(1, np);

% Atmospheric data arrays
atmhPa1             = zeros(1, np);
atmhPa2             = zeros(1, np);
atmwindspeed        = zeros(1, np);
atmwinddirection    = zeros(1, np);
atmrelativehumidity = zeros(1, np);

% Fill data arrays
ctddate = datetime(rsk_processed.profiles.downcast.tstart(profiles_to_read), 'ConvertFrom', 'datenum');
[spcol] = getchannelindex(rsk_processed, {'Sea Pressure'});

if opts.do_binning
    ctdseapres = rsk_processed.data(1).values(:, spcol);
    ctddepth   = gsw_z_from_p(ctdseapres, mean(station_info.lat));
else
    mean_lat = mean(station_info.lat, 'omitnan');
end

[ccol, tcol, scol, vcol, dacol] = getchannelindex(rsk_processed, ...
    {'Conductivity', temp_channel, 'Salinity', 'Velocity', 'Density Anomaly'});

for ip = profiles_to_process
    ctdid(ip) = string(station_info.ID{ip});
    profind = getdataindex(rsk_processed, 'direction', opts.cast, 'profile', profiles_to_read(ip));
    if isfield(rsk_processed.data(profind), 'latitude') && ~isempty(rsk_processed.data(profind).latitude)
        ctdlat(ip) = rsk_processed.data(profind).latitude;
        ctdlon(ip) = rsk_processed.data(profind).longitude;
    else
        ctdlat(ip) = station_info.lat(min(ip, length(station_info.lat)));
        ctdlon(ip) = station_info.lon(min(ip, length(station_info.lon)));
    end

    ns = size(rsk_processed.data(profind).values, 1);
    n_samples(ip) = ns;

    ctdcond(1:ns, ip)    = rsk_processed.data(profind).values(:, ccol);
    ctdtemp(1:ns, ip)    = rsk_processed.data(profind).values(:, tcol);
    ctdsal(1:ns, ip)     = rsk_processed.data(profind).values(:, scol);
    ctdsigma(1:ns, ip)   = rsk_processed.data(profind).values(:, dacol);
    ctdvelprof(1:ns, ip) = rsk_processed.data(profind).values(:, vcol);

    if ~opts.do_binning
        ctdseapres(1:ns, ip)   = rsk_processed.data(profind).values(:, spcol);
        ctdsampltime(1:ns, ip) = rsk_processed.data(profind).tstamp;
    end

    % Soak metadata
    soak_duration_s(ip) = soak_info(ip).duration_s;
    soak_depth_dbar(ip) = soak_info(ip).depth_dbar;
    soak_n_filtered(ip) = soak_info(ip).n_filtered;

    % Atmospheric data (nearest FerryBox hourly mean)
    if ~isempty(atm_info)
        [~, idx] = min(abs(atm_info.data.hourlyMean.DateTime - ctddate(ip)));
        atmhPa1(ip)             = atm_info.data.hourlyMean.AtmPress1_hPa(idx);
        atmhPa2(ip)             = atm_info.data.hourlyMean.AtmPress2_hPa(idx);
        atmwindspeed(ip)        = atm_info.data.hourlyMean.WindSpeed_ms(idx);
        atmwinddirection(ip)    = atm_info.data.hourlyMean.WindDirection(idx);
        atmrelativehumidity(ip) = atm_info.data.hourlyMean.RelativeHumidity(idx);
    end
end

if opts.do_binning
    % ctddepth is 1D (computed above from common pressure grid)
else
    ctddepth = NaN(m, np);
    for ip = profiles_to_process
        ns = n_samples(ip);
        ctddepth(1:ns, ip) = abs(gsw_z_from_p(ctdseapres(1:ns, ip), ctdlat(ip)));
    end
end

% QC flags
ctdseapres_qc = qc_flags.seapres_qc;
ctdtemp_qc    = qc_flags.temp_qc;
ctdcond_qc    = qc_flags.cond_qc;
ctdsal_qc     = qc_flags.sal_qc;
ctdsigma_qc   = qc_flags.sigma_qc;

% Validation decisions (per profile)
validation_decisions = cell(1, np);
for ip = 1:np
    validation_decisions{ip} = validation.profiles(ip).decision;
end

% Output filename: bin025 or native depending on do_binning mode
if opts.do_binning
    outfile = sprintf('%s/PROC_CTD_%s_bin025.mat', opts.output_path, station_info.campaign);
    save(outfile, ...
        'ctddate', 'ctdid', 'ctdlat', 'ctdlon', 'ctddepth', 'ctdseapres', ...
        'ctdcond', 'ctdtemp', 'ctdsal', 'ctdsigma', 'ctdvelprof', ...
        'ctdseapres_qc', 'ctdtemp_qc', 'ctdcond_qc', 'ctdsal_qc', 'ctdsigma_qc', ...
        'n_samples', 'soak_duration_s', 'soak_depth_dbar', 'soak_n_filtered', ...
        'atmhPa1', 'atmhPa2', 'atmwindspeed', 'atmwinddirection', 'atmrelativehumidity', ...
        'validation', 'validation_decisions', 'qc_summary');
else
    outfile = sprintf('%s/PROC_CTD_%s.mat', opts.output_path, station_info.campaign);
    save(outfile, ...
        'ctddate', 'ctdid', 'ctdlat', 'ctdlon', 'ctddepth', 'ctdseapres', 'ctdsampltime', ...
        'ctdcond', 'ctdtemp', 'ctdsal', 'ctdsigma', 'ctdvelprof', ...
        'ctdseapres_qc', 'ctdtemp_qc', 'ctdcond_qc', 'ctdsal_qc', 'ctdsigma_qc', ...
        'n_samples', 'soak_duration_s', 'soak_depth_dbar', 'soak_n_filtered', ...
        'atmhPa1', 'atmhPa2', 'atmwindspeed', 'atmwinddirection', 'atmrelativehumidity', ...
        'validation', 'validation_decisions', 'qc_summary');
end

fprintf('  File saved: %s\n', outfile);

%% 7bis. Soak validation figures (4-panel)
disp(' ')
disp('════════════════════════════════════════════════════════════════════')
disp('  STEP 7bis: SOAK VALIDATION FIGURES (4-panel)')
disp('════════════════════════════════════════════════════════════════════')

if opts.save_figures

    % Colors
    col_soak = [0.85 0.15 0.15];  % red   (L0 soak phase)
    col_l0   = [0.50 0.50 0.50];  % grey  (L0 downcast)
    col_proc = [0.07 0.55 0.32];  % green (processed)
    col_vel  = [0.90 0.45 0.05];  % orange (velocity)
    col_pres = [0.10 0.10 0.10];  % black (pressure)

    sp_zoom = 7;  % dbar — max depth for vertical profile zoom (left panels)
    sp_tmax = 7;  % dbar — time zoom to this depth + 5 s (right panels)

    [spcol_p] = getchannelindex(rsk_processed, 'Sea Pressure');
    [tcol_p]  = getchannelindex(rsk_processed, temp_channel);
    [ccol_p]  = getchannelindex(rsk_processed, 'Conductivity');
    [scol_p]  = getchannelindex(rsk_processed, 'Salinity');
    [vcol_p]  = getchannelindex(rsk_processed, 'Velocity');

    for ip = profiles_to_process
        if iscell(station_info.ID)
            prof_id = station_info.ID{ip};
        else
            prof_id = char(station_info.ID(ip));
        end

        % --- RAW data (reloaded from file for clean baseline) ---
        rsk_raw_v = RSKopen(rsk_file);
        rsk_raw_v = RSKreadprofiles(rsk_raw_v, 'profile', profiles_to_read(ip), 'direction', opts.cast);
        if is_temp1
            raw_ch = {rsk_raw_v.channels.longName};
            idx_t1 = find(strcmpi(raw_ch, 'Temperature1'));
            if ~isempty(idx_t1), rsk_raw_v.channels(idx_t1).longName = 'Temperature'; end
        end
        rsk_raw_v = RSKderiveseapressure(rsk_raw_v, 'patm', opts.patm);
        rsk_raw_v = RSKderivesalinity(rsk_raw_v);
        rsk_raw_v = RSKderivedepth(rsk_raw_v);
        rsk_raw_v = RSKderivevelocity(rsk_raw_v);

        SPcol_rv = getchannelindex(rsk_raw_v, 'Sea Pressure');
        Tcol_rv  = getchannelindex(rsk_raw_v, temp_channel);
        Ccol_rv  = getchannelindex(rsk_raw_v, 'Conductivity');
        Scol_rv  = getchannelindex(rsk_raw_v, 'Salinity');
        Vcol_rv  = getchannelindex(rsk_raw_v, 'Velocity');

        raw_SP = rsk_raw_v.data(1).values(:, SPcol_rv);
        raw_T  = rsk_raw_v.data(1).values(:, Tcol_rv);
        raw_C  = rsk_raw_v.data(1).values(:, Ccol_rv);
        raw_S  = rsk_raw_v.data(1).values(:, Scol_rv);
        raw_V  = rsk_raw_v.data(1).values(:, Vcol_rv);
        raw_t  = rsk_raw_v.data(1).tstamp;
        t_rel  = (raw_t - raw_t(1)) * 86400;  % seconds since profile start

        % --- Processed data ---
        profind_p = getdataindex(rsk_processed, 'direction', opts.cast, 'profile', profiles_to_read(ip));
        proc_SP = rsk_processed.data(profind_p).values(:, spcol_p);
        proc_T  = rsk_processed.data(profind_p).values(:, tcol_p);
        proc_C  = rsk_processed.data(profind_p).values(:, ccol_p);
        proc_S  = rsk_processed.data(profind_p).values(:, scol_p);
        proc_t  = rsk_processed.data(profind_p).tstamp;
        t_rel_p = (proc_t - raw_t(1)) * 86400;  % relative to RAW profile start

        % --- Soak metadata ---
        soak_dur   = soak_info(ip).duration_s;
        soak_dep   = soak_info(ip).depth_dbar;
        soak_nfilt = soak_info(ip).n_filtered;
        soak_pct   = 100 * soak_nfilt / length(raw_SP);

        % --- Split soak / downcast in RAW ---
        n_raw   = length(raw_SP);
        is_soak = false(n_raw, 1);
        is_soak(1:min(soak_nfilt, n_raw)) = true;

        % --- SSS water sample (optional) ---
        % Add station_info.sss = [val1, val2, ...] in proc_run_CTD_by_stations.m
        SSS_val = NaN;
        if isfield(station_info, 'sss') && numel(station_info.sss) >= ip
            SSS_val = station_info.sss(ip);
        end

        % --- Valid data masks ---
        valid_rS = ~isnan(raw_SP) & ~isnan(raw_S);
        valid_rT = ~isnan(raw_SP) & ~isnan(raw_T);
        valid_rC = ~isnan(raw_SP) & ~isnan(raw_C);
        valid_rV = ~isnan(raw_SP) & ~isnan(raw_V);
        valid_pS = ~isnan(proc_SP) & ~isnan(proc_S);
        valid_pT = ~isnan(proc_SP) & ~isnan(proc_T);
        valid_pC = ~isnan(proc_SP) & ~isnan(proc_C);

        % --- X-limits for time panels (full profile) ---
        t_xlim = [0, t_rel(end) + 5];

        % =====================================================================
        % 4-panel figure (2×2)
        % =====================================================================
        fig_v = figure('Position', [50 50 1500 900], 'Color', 'w');

        % ── Top-left: Salinity vs SP, zoom 0–5 dbar ─────────────────────────
        ax1 = subplot(2, 2, 1);
        hold on;

        mk = valid_rS & is_soak & raw_SP <= sp_zoom;
        if any(mk), plot(raw_S(mk), raw_SP(mk), '.', 'Color', col_soak, 'MarkerSize', 12); end
        mk = valid_rS & ~is_soak & raw_SP <= sp_zoom;
        if any(mk), plot(raw_S(mk), raw_SP(mk), '.', 'Color', [0.55 0.70 0.95], 'MarkerSize', 12); end
        mk = valid_pS & proc_SP <= sp_zoom;
        if any(mk), plot(proc_S(mk), proc_SP(mk), '.', 'Color', [0.10 0.30 0.80], 'MarkerSize', 12); end
        if ~isnan(SSS_val)
            plot(SSS_val, 0.2, 'p', 'Color', [0 0 0], ...
                'MarkerFaceColor', [1.0 0.85 0.0], 'MarkerSize', 20, 'LineWidth', 2.0);
        end
        if ~isnan(soak_dep) && soak_dep <= sp_zoom
            yline(soak_dep, '--', 'Color', col_soak, 'LineWidth', 2.0, ...
                'Label', sprintf('%.2f dbar', soak_dep), 'FontSize', 9, ...
                'LabelHorizontalAlignment', 'left');
        end
        hl1 = plot(nan, nan, 'o', 'Color', col_soak, 'MarkerFaceColor', col_soak, 'MarkerSize', 11);
        hl2 = plot(nan, nan, 'o', 'Color', [0.55 0.70 0.95], 'MarkerFaceColor', [0.55 0.70 0.95], 'MarkerSize', 11);
        hl3 = plot(nan, nan, 'o', 'Color', [0.10 0.30 0.80],  'MarkerFaceColor', [0.10 0.30 0.80],  'MarkerSize', 11);
        if ~isnan(SSS_val)
            hl4 = plot(nan, nan, 'p', 'Color', [0 0 0], 'MarkerFaceColor', [1 0.85 0], 'MarkerSize', 14);
            legend([hl1 hl2 hl3 hl4], {'SOAK', 'RAW', 'PROC', 'AutoSal (0.2 m)'}, ...
                'Location', 'best', 'FontSize', 9);
        else
            legend([hl1 hl2 hl3], {'SOAK', 'RAW', 'PROC'}, ...
                'Location', 'best', 'FontSize', 9);
        end
        set(ax1, 'YDir', 'reverse', 'FontSize', 11);
        ylim([-0.5 sp_zoom]);
        xlabel('Practical Salinity (PSU)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Sea Pressure (dbar)', 'FontSize', 11, 'FontWeight', 'bold');
        title('Salinity – zoom 0-7 dbar', 'FontSize', 10);
        grid on; box on;

        % ── Bottom-left: Temperature vs SP, zoom 0–5 dbar ───────────────────
        ax2 = subplot(2, 2, 3);
        hold on;

        mk = valid_rT & is_soak & raw_SP <= sp_zoom;
        if any(mk), plot(raw_T(mk), raw_SP(mk), '.', 'Color', col_soak, 'MarkerSize', 12); end
        ht_s = plot(nan, nan, 'o', 'Color', col_soak,        'MarkerFaceColor', col_soak,        'MarkerSize', 11);
        mk = valid_rT & is_soak & raw_SP <= sp_zoom;
        if any(mk), plot(raw_T(mk), raw_SP(mk), '.', 'Color', col_soak, 'MarkerSize', 12); end
        ht_r = plot(nan, nan, 'o', 'Color', [0.55 0.85 0.55], 'MarkerFaceColor', [0.55 0.85 0.55], 'MarkerSize', 11);
        mk = valid_rT & ~is_soak & raw_SP <= sp_zoom;
        if any(mk), plot(raw_T(mk), raw_SP(mk), '.', 'Color', [0.55 0.85 0.55], 'MarkerSize', 12); end
        ht_p = plot(nan, nan, 'o', 'Color', col_proc,         'MarkerFaceColor', col_proc,         'MarkerSize', 11);
        mk = valid_pT & proc_SP <= sp_zoom;
        if any(mk), plot(proc_T(mk), proc_SP(mk), '.', 'Color', col_proc, 'MarkerSize', 12); end
        if ~isnan(soak_dep) && soak_dep <= sp_zoom
            yline(soak_dep, '--', 'Color', col_soak, 'LineWidth', 2.0, ...
                'Label', sprintf('%.2f dbar', soak_dep), 'FontSize', 9, ...
                'LabelHorizontalAlignment', 'left');
        end
        legend([ht_s ht_r ht_p], {'SOAK', 'RAW', 'PROC'}, 'Location', 'best', 'FontSize', 9);
        set(ax2, 'YDir', 'reverse', 'FontSize', 11);
        ylim([-0.5 sp_zoom]);
        xlabel('Temperature (°C)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Sea Pressure (dbar)', 'FontSize', 11, 'FontWeight', 'bold');
        title('Temperature – zoom 0-7 dbar', 'FontSize', 10);
        grid on; box on;

        % ── Top-right: Temperature and Salinity vs time ──────────────────────
        ax3 = subplot(2, 2, 2);
        hold on;

        yyaxis left
        ht1 = plot(t_rel(valid_rT),   raw_T(valid_rT),   '.', ...
            'Color', [0.55 0.85 0.55], 'MarkerSize', 8);
        ht2 = plot(t_rel_p(valid_pT), proc_T(valid_pT),  '.', ...
            'Color', col_proc, 'MarkerSize', 8);
        ylabel('Temperature (°C)', 'FontSize', 11, 'FontWeight', 'bold');
        set(ax3.YAxis(1), 'Color', col_proc);

        yyaxis right
        hs1 = plot(t_rel(valid_rS),   raw_S(valid_rS),   '.', ...
            'Color', [0.55 0.70 0.95], 'MarkerSize', 8);
        hs2 = plot(t_rel_p(valid_pS), proc_S(valid_pS),  '.', ...
            'Color', [0.10 0.30 0.80], 'MarkerSize', 8);
        ylabel('Salinity (PSU)', 'FontSize', 11, 'FontWeight', 'bold');
        set(ax3.YAxis(2), 'Color', [0.10 0.30 0.80]);

        if ~isnan(soak_dur) && soak_dur <= t_xlim(2)
            xline(soak_dur, '-', sprintf('Soak: %.0f s', soak_dur), ...
                'Color', col_soak, 'LineWidth', 2.5, 'FontSize', 9, ...
                'LabelVerticalAlignment', 'bottom');
        end
        legend([ht1 ht2 hs1 hs2], {'T RAW', 'T PROC', 'S RAW', 'S PROC'}, ...
            'Location', 'best', 'FontSize', 9);
        xlim(t_xlim);
        xlabel('Time since profile start (s)', 'FontSize', 11, 'FontWeight', 'bold');
        title('Temperature & Salinity vs Time', 'FontSize', 10);
        grid on; box on;

        % ── Bottom-right: Velocity + Sea Pressure vs time ────────────────────
        ax4 = subplot(2, 2, 4);
        hold on;

        yyaxis left
        soak_x_end = min(soak_dur, t_xlim(2));
        if ~isnan(soak_dur) && soak_dur > 0
            patch([0 soak_x_end soak_x_end 0], [-0.3 -0.3 1.5 1.5], ...
                col_soak, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
        end
        h_vel = plot(t_rel(valid_rV), raw_V(valid_rV), '-', ...
            'Color', col_vel, 'LineWidth', 2.0);
        yline(opts.velocity_threshold_loop, '--', 'Color', col_vel, 'LineWidth', 1.5, ...
            'Label', sprintf('%.2f m/s (loop)', opts.velocity_threshold_loop), 'FontSize', 9, ...
            'LabelVerticalAlignment', 'bottom');
        ylim([-0.3 1.5]);
        ylabel('Descent velocity (m/s)', 'FontSize', 11, 'FontWeight', 'bold');
        set(ax4.YAxis(1), 'Color', col_vel);

        yyaxis right
        h_pres = plot(t_rel(~isnan(raw_SP)), raw_SP(~isnan(raw_SP)), '-', ...
            'Color', col_pres, 'LineWidth', 2.0);
        set(gca, 'YDir', 'reverse');
        ylim([-0.5 max(raw_SP(~isnan(raw_SP))) * 1.05]);
        ylabel('Sea Pressure (dbar)', 'FontSize', 11, 'FontWeight', 'bold');
        set(ax4.YAxis(2), 'Color', col_pres);

        if ~isnan(soak_dur) && soak_dur > 0 && soak_dur <= t_xlim(2)
            xline(soak_dur, '-', sprintf('Soak end: %.0f s  (%.2f dbar)', soak_dur, soak_dep), ...
                'Color', col_soak, 'LineWidth', 2.5, 'FontSize', 9, ...
                'LabelVerticalAlignment', 'bottom');
        end
        xlim(t_xlim);
        xlabel('Time since profile start (s)', 'FontSize', 11, 'FontWeight', 'bold');
        title('Velocity and Sea Pressure vs Time', 'FontSize', 10);
        legend([h_vel h_pres], {'Velocity', 'Sea Pressure'}, 'Location', 'best', 'FontSize', 9);
        box on;

        % --- Overall title ---
        freq_str = '16 Hz';
        if ~is_highfreq, freq_str = '2 Hz'; end
        sgtitle(sprintf('L0 vs Processed — %s  (%s, %s)  |  Soak: %.0f s @ %.2f dbar  (%.1f%%)', ...
            prof_id, station_info.campaign, freq_str, soak_dur, soak_dep, soak_pct), ...
            'FontSize', 12, 'FontWeight', 'bold');

        % --- Save figure ---
        figname_sv = sprintf('%s/SoakValidation_%s_profile%02d_%s.png', ...
                            opts.fig_path, station_info.campaign, ip, prof_id);
        exportgraphics(fig_v, figname_sv, 'Resolution', 200, 'BackgroundColor', 'white');
        close(fig_v);
        fprintf('  Soak validation: profile %02d (%s)\n', ip, prof_id);
    end
    disp('  [OK] Soak validation figures generated')
else
    disp('  [SKIP] save_figures = false')
end

%% 8. Final summary
disp(' ')
disp('════════════════════════════════════════════════════════════════════')
disp('  SUMMARY — PROCESSING COMPLETE')
disp('════════════════════════════════════════════════════════════════════')
disp(' ')
fprintf('  Campaign         : %s\n', station_info.campaign);
fprintf('  Total profiles   : %d\n', np);
fprintf('  Processed        : %d\n', np_process);
fprintf('  Excluded (SKIP)  : %d\n', length(validation.skipped));
fprintf('  Suspect profiles : %d\n', length(validation.suspects));
disp(' ')
disp('  Integrated validation:')
n_ok      = sum(strcmp(validation_decisions, 'OK'));
n_suspect = sum(strcmp(validation_decisions, 'SUSPECT'));
n_bad     = sum(strcmp(validation_decisions, 'BAD'));
n_skip    = sum(strcmp(validation_decisions, 'SKIP'));
fprintf('    OK      : %d\n', n_ok);
fprintf('    SUSPECT : %d\n', n_suspect);
fprintf('    BAD     : %d\n', n_bad);
fprintf('    SKIP    : %d\n', n_skip);
disp(' ')
disp('  Quality Control (QC):')
fprintf('    Tests applied : %d\n',    length(qc_summary.tests_applied));
fprintf('    Temp — Good   : %.1f%%\n', qc_summary.flags_count.temp.pct_good);
fprintf('    Sal  — Good   : %.1f%%\n', qc_summary.flags_count.sal.pct_good);
disp(' ')
fprintf('  Output file : %s\n', outfile);
if opts.save_figures
    fprintf('  Figures     : %s\n', opts.fig_path);
end
disp('════════════════════════════════════════════════════════════════════')

%% Output struct
output.matfile      = outfile;
output.np_total     = np;
output.np_processed = np_process;
output.np_skipped   = length(validation.skipped);
output.validation   = validation;
output.qc_summary   = qc_summary;

end
