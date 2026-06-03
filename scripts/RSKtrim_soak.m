function [rsk, soak_info] = RSKtrim_soak_auto(rsk, varargin)
%% RSKtrim_soak_auto - Automatic soak phase detection and removal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nina HOAREAU, January 2025 — updated May 2026 (BPL, ICM-CSIC)
% ---
% Automatically detects the soak phase (stagnation at surface) at the start
% of each CTD profile by analysing the descent velocity, then removes those
% data points.
%
% TWO EXCLUSIVE MODES (selected via 'method'):
%
%   MODE 'fixed_time' (DEFAULT):
%     Removes data from the start of the profile up to fixed_time_s seconds.
%     Simple, reproducible, recommended.
%     Default: fixed_time_s = 30 seconds.
%
%   MODE 'velocity':
%     Automatic detection based on descent velocity:
%     1. Skip the first 5 seconds (physical entry of the sensor into the water).
%        At 16 Hz: skip = 81 pts → t_rel(81) = 5.0 s. At 2 Hz: skip = 11 pts.
%     2. From that point, find the first window of min_consecutive consecutive
%        points where V_smooth > velocity_threshold.
%        This marks the start of the true descent → soak_end_idx.
%     If no window found → soak_end_idx = skip_pts (minimum 5 s removed).
%
% NOTE: vel_static is accepted as a parameter for backward compatibility with
%   older calling scripts, but is not used in this algorithm.
%
% USAGE:
%   % Fixed-time mode (default, 30 s):
%   [rsk, soak_info] = RSKtrim_soak_auto(rsk)
%   [rsk, soak_info] = RSKtrim_soak_auto(rsk, 'fixed_time_s', 20)
%
%   % Velocity mode:
%   [rsk, soak_info] = RSKtrim_soak_auto(rsk, 'method', 'velocity')
%   [rsk, soak_info] = RSKtrim_soak_auto(rsk, 'method', 'velocity', ...
%                          'velocity_threshold', 0.15)
%
% INPUTS:
%   rsk - RSK structure after RSKreadprofiles, RSKderiveseapressure,
%         RSKderivedepth and RSKderivevelocity
%
% OPTIONAL PARAMETERS:
%   'method'             - Detection mode: 'fixed_time' (default) or 'velocity'
%   'fixed_time_s'       - [fixed_time mode] Fixed soak duration in seconds.
%                          Default: 30 s
%   'velocity_threshold' - [velocity mode] Speed threshold for "descending"
%                          Default: 0.25 m/s
%   'vel_static'         - Accepted for backward compatibility, not used.
%                          Default: 0.10 m/s
%   'window_size'        - [velocity mode] Window for moving average smoothing.
%                          Default: 16 points (~1 s at 16 Hz)
%   'min_consecutive'    - [velocity mode] Consecutive points required to confirm
%                          descent. Default: 32 points (~2 s at 16 Hz)
%   'action'             - 'nan' (replace with NaN) or 'remove' (delete rows).
%                          Default: 'nan'
%   'min_soak_depth'     - Minimum soak depth (Sea Pressure, dbar).
%                          If detected soak end < min_soak_depth, extend to
%                          that depth. Default: 0.5 dbar
%   'verbose'            - Print summary information (true/false).
%                          Default: true
%
% OUTPUTS:
%   rsk       - RSK structure with soak data removed/NaN-ed
%   soak_info - Structure with soak information for each profile:
%               .profile     - Profile number
%               .duration_s  - Soak duration (seconds)
%               .depth_dbar  - Sea pressure at soak end (dbar)
%               .n_filtered  - Number of points removed/NaN-ed
%               .method      - Method used ('fixed_time' or 'velocity')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parse inputs
p = inputParser;
addRequired(p, 'rsk', @isstruct);
addParameter(p, 'method',             'fixed_time', @ischar);
addParameter(p, 'fixed_time_s',              30,    @isnumeric);
addParameter(p, 'velocity_threshold',      0.25,    @isnumeric);
addParameter(p, 'vel_static',              0.10,    @isnumeric);  % backward compat, not used
addParameter(p, 'window_size',               16,    @isnumeric);
addParameter(p, 'min_consecutive',           32,    @isnumeric);
addParameter(p, 'action',               'nan',      @ischar);
addParameter(p, 'min_soak_depth',           0.5,    @isnumeric);
addParameter(p, 'verbose',              true,       @islogical);

parse(p, rsk, varargin{:});

method             = p.Results.method;
fixed_time_s       = p.Results.fixed_time_s;
velocity_threshold = p.Results.velocity_threshold;
window_size        = p.Results.window_size;
min_consecutive    = p.Results.min_consecutive;
action             = p.Results.action;
min_soak_depth     = p.Results.min_soak_depth;
verbose            = p.Results.verbose;

if ~ismember(method, {'fixed_time', 'velocity'})
    error('RSKtrim_soak_auto: ''method'' must be ''fixed_time'' or ''velocity''. Got: ''%s''.', method);
end

%% Check required channels
try
    SPcol = getchannelindex(rsk, 'Sea Pressure');
catch
    error('RSKtrim_soak_auto: Sea Pressure channel not found. Run RSKderiveseapressure first.');
end

if strcmp(method, 'velocity')
    try
        Vcol = getchannelindex(rsk, 'Velocity');
    catch
        error('RSKtrim_soak_auto: Velocity channel not found. Run RSKderivedepth and RSKderivevelocity first.');
    end
else
    Vcol = [];  % not used in fixed_time mode
end

%% Number of profiles
np = length(rsk.data);

if verbose
    disp('─────────────────────────────────────────────────────')
    disp('  RSKtrim_soak_auto - Automatic soak detection')
    disp('─────────────────────────────────────────────────────')
    fprintf('  Mode : %s\n', upper(method))
    if strcmp(method, 'fixed_time')
        fprintf('    fixed_time_s         : %.1f s\n', fixed_time_s)
    else
        fprintf('    velocity_threshold   : %.2f m/s\n', velocity_threshold)
        fprintf('    window_size          : %d points\n', window_size)
        fprintf('    min_consecutive      : %d points\n', min_consecutive)
    end
    fprintf('    min_soak_depth       : %.2f dbar\n', min_soak_depth)
    fprintf('    action               : %s\n', action)
    disp(' ')
    fprintf('  %-8s %10s %10s %10s\n', 'Profile', 'Duration(s)', 'SP_soak', 'Pts removed');
    disp('  ────────────────────────────────────────────')
end

%% Initialise output structure
soak_info = struct('profile', {}, 'duration_s', {}, 'depth_dbar', {}, ...
                   'n_filtered', {}, 'method', {});

%% Loop over profiles
for ip = 1:np

    SP   = rsk.data(ip).values(:, SPcol);
    time = rsk.data(ip).tstamp;

    n_pts = length(SP);

    % Relative time (seconds from profile start)
    t_rel = (time - time(1)) * 24 * 60 * 60;

    % ─────────────────────────────────────────────────────────────────────
    % FIXED-TIME MODE
    % ─────────────────────────────────────────────────────────────────────
    if strcmp(method, 'fixed_time')

        idx_fixed = find(t_rel >= fixed_time_s, 1, 'first');

        if isempty(idx_fixed)
            % fixed_time_s exceeds total profile duration — safety fallback
            idx_fixed = n_pts;
            warning('RSKtrim_soak_auto: profile %d — fixed_time_s (%.1f s) exceeds total profile duration. Entire profile removed.', ip, fixed_time_s);
        end

        soak_end_idx = idx_fixed;

    % ─────────────────────────────────────────────────────────────────────
    % VELOCITY MODE
    % ─────────────────────────────────────────────────────────────────────
    else

        V        = rsk.data(ip).values(:, Vcol);
        V_smooth = movmean(V, window_size, 'omitnan');

        dt           = (time(2) - time(1)) * 86400;        % seconds between samples
        skip_pts     = max(1, floor(5.0 / dt) + 1);        % 81 pts at 16 Hz, 11 pts at 2 Hz
        search_start = min(skip_pts, n_pts - min_consecutive);

        soak_end_idx = search_start;   % default: minimum 5 s removed

        for i = search_start:(n_pts - min_consecutive)
            if all(V_smooth(i:i+min_consecutive-1) > velocity_threshold)
                soak_end_idx = i;
                break;
            end
        end

    end

    % ─────────────────────────────────────────────────────────────────────
    % Apply minimum soak depth (common to both modes)
    % ─────────────────────────────────────────────────────────────────────

    if SP(soak_end_idx) < min_soak_depth
        idx_min_depth = find(SP >= min_soak_depth, 1, 'first');
        if ~isempty(idx_min_depth) && idx_min_depth > soak_end_idx
            soak_end_idx = idx_min_depth;
        end
    end

    % ─────────────────────────────────────────────────────────────────────
    % Apply filter
    % ─────────────────────────────────────────────────────────────────────

    nan_until = soak_end_idx;

    if strcmp(action, 'nan')
        rsk.data(ip).values(1:nan_until, :) = NaN;
    elseif strcmp(action, 'remove')
        rsk.data(ip).values(1:nan_until, :) = [];
        rsk.data(ip).tstamp(1:nan_until)    = [];
    end

    % ─────────────────────────────────────────────────────────────────────
    % Store soak info
    % ─────────────────────────────────────────────────────────────────────

    soak_duration = t_rel(nan_until);
    soak_depth    = SP(nan_until);
    n_filtered    = nan_until;

    soak_info(ip).profile    = ip;
    soak_info(ip).duration_s = soak_duration;
    soak_info(ip).depth_dbar = soak_depth;
    soak_info(ip).n_filtered = n_filtered;
    soak_info(ip).method     = method;

    if verbose
        fprintf('  %-8d %10.1f %10.2f %10d\n', ip, soak_duration, soak_depth, n_filtered);
    end

end

%% Summary
if verbose
    disp('  ────────────────────────────────────────────')
    disp(' ')
    disp('  Summary:')
    durations = [soak_info.duration_s];
    depths    = [soak_info.depth_dbar];
    fprintf('    Mean soak duration   : %.1f s (%.1f - %.1f)\n', ...
            mean(durations), min(durations), max(durations));
    fprintf('    Mean soak depth      : %.2f dbar (%.2f - %.2f)\n', ...
            mean(depths), min(depths), max(depths));
    disp(' ')
    disp('─────────────────────────────────────────────────────')
end

end
