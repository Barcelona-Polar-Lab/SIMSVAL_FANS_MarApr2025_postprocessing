function [qc_flags, qc_summary] = apply_QC_tests(data, varargin)
%% APPLY_QC_TESTS - Apply Quality Control tests to CTD data (post-binning)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nina HOAREAU, January 2025 (BPL, ICM-CSIC)
% ---
% This function applies a series of QC tests to binned CTD data following
% international standards (QARTOD, Argo, SeaDataNet).
%
% QC FLAGS (SeaDataNet convention):
%   0 = No QC performed
%   1 = Good data
%   2 = Probably good data
%   3 = Suspect/Doubtful data
%   4 = Bad data
%   5 = Changed value
%   8 = Interpolated value
%   9 = Missing data
%
% USAGE:
%   [qc_flags, qc_summary] = apply_QC_tests(data)
%   [qc_flags, qc_summary] = apply_QC_tests(data, 'region', 'arctic')
%   [qc_flags, qc_summary] = apply_QC_tests(data, 'param', value, ...)
%
% INPUT:
%   data : structure with fields:
%       .seapres   (m x np) - Sea pressure [dbar]
%       .temp      (m x np) - Temperature [°C]
%       .cond      (m x np) - Conductivity [mS/cm]
%       .sal       (m x np) - Practical Salinity [PSU]
%       .sigma     (m x np) - Density anomaly [kg/m³]
%       .lat       (1 x np) - Latitude [°N]
%       .lon       (1 x np) - Longitude [°E]
%
% OUTPUT:
%   qc_flags : structure with QC flag matrices for each variable
%       .seapres_qc  (m x np)
%       .temp_qc     (m x np)
%       .cond_qc     (m x np)
%       .sal_qc      (m x np)
%       .sigma_qc    (m x np)
%
%   qc_summary : structure with summary statistics per profile
%
% OPTIONAL PARAMETERS (name-value pairs):
%   'region'          - 'global' or 'arctic' (default: 'arctic')
%   'flatline_n'      - Number of consecutive identical values for flat line (default: 5)
%   'density_tol'     - Tolerance for density inversion [kg/m³] (default: 0.03)
%   'verbose'         - Display progress messages (default: true)
%
% TESTS IMPLEMENTED:
%   1. NaN Check
%   2. Gross Range Test (global limits)
%   3. Regional Range Test (Arctic-specific)
%   4. Flat Line Test
%   5. Gradient Test (vertical)
%   6. Density Inversion Test
%   7. Pressure Monotonicity Check
%
% REFERENCES:
%   - QARTOD Manual for T&S (IOOS, 2020)
%   - Argo Quality Control Manual (2023)
%   - SeaDataNet Data Quality Control Procedures V2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parse input arguments
p = inputParser;
addRequired(p, 'data', @isstruct);
addParameter(p, 'region', 'arctic', @ischar);
addParameter(p, 'flatline_n', 5, @isnumeric);
addParameter(p, 'density_tol', 0.03, @isnumeric);
addParameter(p, 'skip_gradient', false, @islogical);   % skip Test 5 for native-resolution data
addParameter(p, 'verbose', true, @islogical);

parse(p, data, varargin{:});
opts = p.Results;

%% Define QC thresholds
% =========================================================================
% GLOBAL THRESHOLDS (from QARTOD/Argo)
% =========================================================================
thresholds.global.temp_min = -2.5;      % °C
thresholds.global.temp_max = 40.0;      % °C
thresholds.global.sal_min = 2.0;        % PSU
thresholds.global.sal_max = 42.0;       % PSU
thresholds.global.cond_min = 0.0;       % mS/cm
thresholds.global.cond_max = 70.0;      % mS/cm
thresholds.global.pres_min = -0.5;      % dbar (small tolerance for surface noise)
thresholds.global.pres_max = 6500.0;    % dbar

% =========================================================================
% ARCTIC THRESHOLDS (regional, more restrictive)
% =========================================================================
thresholds.arctic.temp_min = -2.0;      % °C (freezing point ~-1.9°C)
thresholds.arctic.temp_max = 15.0;      % °C
thresholds.arctic.sal_min = 20.0;       % PSU (fjords can have fresher water)
thresholds.arctic.sal_max = 36.0;       % PSU
thresholds.arctic.cond_min = 5.0;       % mS/cm
thresholds.arctic.cond_max = 50.0;      % mS/cm
thresholds.arctic.pres_min = 0.0;       % dbar
thresholds.arctic.pres_max = 500.0;     % dbar (shallow Arctic shelf/fjords)

% =========================================================================
% GRADIENT THRESHOLDS (vertical, per dbar)
% =========================================================================
thresholds.grad.temp = 2.0;             % °C/dbar (suspect if exceeded)
thresholds.grad.sal = 1.0;              % PSU/dbar
thresholds.grad.cond = 5.0;             % mS/cm/dbar

% =========================================================================
% OTHER PARAMETERS
% =========================================================================
thresholds.flatline_n = opts.flatline_n;
thresholds.density_tol = opts.density_tol;

%% Get data dimensions
[m, np] = size(data.temp);

if opts.verbose
    fprintf('\n========================================\n');
    fprintf('  APPLY_QC_TESTS - CTD Quality Control\n');
    fprintf('========================================\n');
    fprintf('  Region: %s\n', upper(opts.region));
    fprintf('  Profiles: %d\n', np);
    fprintf('  Depth levels: %d\n', m);
    fprintf('========================================\n\n');
end

%% Initialize QC flag matrices (0 = no QC performed)
qc_flags.seapres_qc = zeros(m, np, 'int8');
qc_flags.temp_qc = zeros(m, np, 'int8');
qc_flags.cond_qc = zeros(m, np, 'int8');
qc_flags.sal_qc = zeros(m, np, 'int8');
qc_flags.sigma_qc = zeros(m, np, 'int8');

% Initialize summary structure
qc_summary.n_profiles = np;
qc_summary.n_levels = m;
qc_summary.tests_applied = {};
qc_summary.flags_count = struct();

%% ========================================================================
% TEST 1: NaN CHECK
% =========================================================================
if opts.verbose, fprintf('TEST 1: NaN Check... '); end

% Flag 9 for missing data (NaN)
qc_flags.seapres_qc(isnan(data.seapres)) = 9;
qc_flags.temp_qc(isnan(data.temp)) = 9;
qc_flags.cond_qc(isnan(data.cond)) = 9;
qc_flags.sal_qc(isnan(data.sal)) = 9;
qc_flags.sigma_qc(isnan(data.sigma)) = 9;

n_nan = sum(isnan(data.temp(:)));
if opts.verbose, fprintf('Done. (%d NaN values flagged)\n', n_nan); end
qc_summary.tests_applied{end+1} = 'NaN Check';

%% ========================================================================
% TEST 2: GROSS RANGE TEST (Global limits)
% =========================================================================
if opts.verbose, fprintf('TEST 2: Gross Range Test (Global)... '); end

% Pressure
idx_bad = data.seapres < thresholds.global.pres_min | ...
          data.seapres > thresholds.global.pres_max;
qc_flags.seapres_qc(idx_bad & qc_flags.seapres_qc == 0) = 4;

% Temperature
idx_bad = data.temp < thresholds.global.temp_min | ...
          data.temp > thresholds.global.temp_max;
qc_flags.temp_qc(idx_bad & qc_flags.temp_qc == 0) = 4;

% Conductivity
idx_bad = data.cond < thresholds.global.cond_min | ...
          data.cond > thresholds.global.cond_max;
qc_flags.cond_qc(idx_bad & qc_flags.cond_qc == 0) = 4;

% Salinity
idx_bad = data.sal < thresholds.global.sal_min | ...
          data.sal > thresholds.global.sal_max;
qc_flags.sal_qc(idx_bad & qc_flags.sal_qc == 0) = 4;

if opts.verbose, fprintf('Done.\n'); end
qc_summary.tests_applied{end+1} = 'Gross Range (Global)';

%% ========================================================================
% TEST 3: REGIONAL RANGE TEST (Arctic-specific)
% =========================================================================
if strcmpi(opts.region, 'arctic')
    if opts.verbose, fprintf('TEST 3: Regional Range Test (Arctic)... '); end

    % Pressure (regional max)
    idx_suspect = data.seapres > thresholds.arctic.pres_max;
    qc_flags.seapres_qc(idx_suspect & qc_flags.seapres_qc == 0) = 3;

    % Temperature
    idx_suspect = data.temp < thresholds.arctic.temp_min | ...
                  data.temp > thresholds.arctic.temp_max;
    qc_flags.temp_qc(idx_suspect & qc_flags.temp_qc == 0) = 3;

    % Conductivity
    idx_suspect = data.cond < thresholds.arctic.cond_min | ...
                  data.cond > thresholds.arctic.cond_max;
    qc_flags.cond_qc(idx_suspect & qc_flags.cond_qc == 0) = 3;

    % Salinity
    idx_suspect = data.sal < thresholds.arctic.sal_min | ...
                  data.sal > thresholds.arctic.sal_max;
    qc_flags.sal_qc(idx_suspect & qc_flags.sal_qc == 0) = 3;

    if opts.verbose, fprintf('Done.\n'); end
    qc_summary.tests_applied{end+1} = 'Regional Range (Arctic)';
end

%% ========================================================================
% TEST 4: FLAT LINE TEST
% =========================================================================
if opts.verbose, fprintf('TEST 4: Flat Line Test (n=%d)... ', thresholds.flatline_n); end

for ip = 1:np
    % Pressure
    qc_flags.seapres_qc(:,ip) = apply_flatline_test(...
        data.seapres(:,ip), qc_flags.seapres_qc(:,ip), thresholds.flatline_n);

    % Temperature
    qc_flags.temp_qc(:,ip) = apply_flatline_test(...
        data.temp(:,ip), qc_flags.temp_qc(:,ip), thresholds.flatline_n);

    % Conductivity
    qc_flags.cond_qc(:,ip) = apply_flatline_test(...
        data.cond(:,ip), qc_flags.cond_qc(:,ip), thresholds.flatline_n);

    % Salinity
    qc_flags.sal_qc(:,ip) = apply_flatline_test(...
        data.sal(:,ip), qc_flags.sal_qc(:,ip), thresholds.flatline_n);
end

if opts.verbose, fprintf('Done.\n'); end
qc_summary.tests_applied{end+1} = sprintf('Flat Line (n=%d)', thresholds.flatline_n);

%% ========================================================================
% TEST 5: GRADIENT TEST (Vertical)
% Note: this test is designed for vertically binned data (e.g. 0.25 dbar bins).
% At native sensor resolution (16 Hz), point-to-point pressure differences
% are ~0.003 dbar, making dT/dP highly sensitive to measurement noise.
% Pass 'skip_gradient', true when calling with native-resolution data.
% =========================================================================
if ~opts.skip_gradient
    if opts.verbose, fprintf('TEST 5: Gradient Test (Vertical)... '); end

    for ip = 1:np
        % Get pressure differences for gradient calculation
        dp = diff(data.seapres(:,ip));
        dp(dp == 0) = NaN;  % Avoid division by zero

        % Temperature gradient
        dT = abs(diff(data.temp(:,ip))) ./ abs(dp);
        idx_suspect = find(dT > thresholds.grad.temp);
        for ii = idx_suspect'
            if qc_flags.temp_qc(ii,ip) == 0, qc_flags.temp_qc(ii,ip) = 3; end
            if qc_flags.temp_qc(ii+1,ip) == 0, qc_flags.temp_qc(ii+1,ip) = 3; end
        end

        % Salinity gradient
        dS = abs(diff(data.sal(:,ip))) ./ abs(dp);
        idx_suspect = find(dS > thresholds.grad.sal);
        for ii = idx_suspect'
            if qc_flags.sal_qc(ii,ip) == 0, qc_flags.sal_qc(ii,ip) = 3; end
            if qc_flags.sal_qc(ii+1,ip) == 0, qc_flags.sal_qc(ii+1,ip) = 3; end
        end

        % Conductivity gradient
        dC = abs(diff(data.cond(:,ip))) ./ abs(dp);
        idx_suspect = find(dC > thresholds.grad.cond);
        for ii = idx_suspect'
            if qc_flags.cond_qc(ii,ip) == 0, qc_flags.cond_qc(ii,ip) = 3; end
            if qc_flags.cond_qc(ii+1,ip) == 0, qc_flags.cond_qc(ii+1,ip) = 3; end
        end
    end

    if opts.verbose, fprintf('Done.\n'); end
    qc_summary.tests_applied{end+1} = 'Gradient (Vertical)';
else
    if opts.verbose, fprintf('TEST 5: Gradient Test — SKIPPED (native-resolution data, not applicable).\n'); end
    qc_summary.tests_applied{end+1} = 'Gradient (Vertical) — SKIPPED (native resolution)';
end

%% ========================================================================
% TEST 6: DENSITY INVERSION TEST
% =========================================================================
if opts.verbose, fprintf('TEST 6: Density Inversion Test (tol=%.3f kg/m³)... ', thresholds.density_tol); end

for ip = 1:np
    dsigma = diff(data.sigma(:,ip));

    % Density should increase with depth (dsigma > 0)
    % Flag as suspect if inversion exceeds tolerance
    idx_inversion = find(dsigma < -thresholds.density_tol);

    for ii = idx_inversion'
        % Flag both points involved in the inversion
        if qc_flags.sigma_qc(ii,ip) == 0, qc_flags.sigma_qc(ii,ip) = 3; end
        if qc_flags.sigma_qc(ii+1,ip) == 0, qc_flags.sigma_qc(ii+1,ip) = 3; end

        % Also flag associated T and S as they contribute to density
        if qc_flags.temp_qc(ii,ip) == 0, qc_flags.temp_qc(ii,ip) = 3; end
        if qc_flags.temp_qc(ii+1,ip) == 0, qc_flags.temp_qc(ii+1,ip) = 3; end
        if qc_flags.sal_qc(ii,ip) == 0, qc_flags.sal_qc(ii,ip) = 3; end
        if qc_flags.sal_qc(ii+1,ip) == 0, qc_flags.sal_qc(ii+1,ip) = 3; end
    end
end

if opts.verbose, fprintf('Done.\n'); end
qc_summary.tests_applied{end+1} = sprintf('Density Inversion (tol=%.3f)', thresholds.density_tol);

%% ========================================================================
% TEST 7: PRESSURE MONOTONICITY CHECK (for downcast data)
% =========================================================================
if opts.verbose, fprintf('TEST 7: Pressure Monotonicity Check... '); end

for ip = 1:np
    dp = diff(data.seapres(:,ip));

    % For downcast, pressure should increase (dp > 0)
    % After binning, this should already be OK, but check anyway
    idx_nonmono = find(dp < 0);

    for ii = idx_nonmono'
        if qc_flags.seapres_qc(ii,ip) == 0, qc_flags.seapres_qc(ii,ip) = 3; end
        if qc_flags.seapres_qc(ii+1,ip) == 0, qc_flags.seapres_qc(ii+1,ip) = 3; end
    end
end

if opts.verbose, fprintf('Done.\n'); end
qc_summary.tests_applied{end+1} = 'Pressure Monotonicity';

%% ========================================================================
% FINAL STEP: Mark remaining zeros as "Good" (flag = 1)
% =========================================================================
if opts.verbose, fprintf('Finalizing: Marking good data... '); end

qc_flags.seapres_qc(qc_flags.seapres_qc == 0) = 1;
qc_flags.temp_qc(qc_flags.temp_qc == 0) = 1;
qc_flags.cond_qc(qc_flags.cond_qc == 0) = 1;
qc_flags.sal_qc(qc_flags.sal_qc == 0) = 1;
qc_flags.sigma_qc(qc_flags.sigma_qc == 0) = 1;

if opts.verbose, fprintf('Done.\n'); end

%% ========================================================================
% COMPUTE SUMMARY STATISTICS
% =========================================================================
qc_summary.flags_count.seapres = count_flags(qc_flags.seapres_qc);
qc_summary.flags_count.temp = count_flags(qc_flags.temp_qc);
qc_summary.flags_count.cond = count_flags(qc_flags.cond_qc);
qc_summary.flags_count.sal = count_flags(qc_flags.sal_qc);
qc_summary.flags_count.sigma = count_flags(qc_flags.sigma_qc);

if opts.verbose
    fprintf('\n========================================\n');
    fprintf('  QC SUMMARY\n');
    fprintf('========================================\n');
    print_qc_summary(qc_summary);
    fprintf('========================================\n\n');
end

end

%% ========================================================================
% HELPER FUNCTIONS
% =========================================================================

function qc = apply_flatline_test(data_vec, qc_vec, n_flat)
% Apply flat line test to a single profile vector
% Flag as suspect (3) if more than n_flat consecutive identical values

    qc = qc_vec;
    n = length(data_vec);

    if n < n_flat
        return;
    end

    count = 1;
    for i = 2:n
        if ~isnan(data_vec(i)) && data_vec(i) == data_vec(i-1)
            count = count + 1;
        else
            if count >= n_flat
                % Flag all values in the flat segment
                for j = (i-count):(i-1)
                    if qc(j) == 0
                        qc(j) = 3;
                    end
                end
            end
            count = 1;
        end
    end

    % Check last segment
    if count >= n_flat
        for j = (n-count+1):n
            if qc(j) == 0
                qc(j) = 3;
            end
        end
    end
end

function counts = count_flags(qc_matrix)
% Count occurrences of each flag value
    counts.flag_1_good = sum(qc_matrix(:) == 1);
    counts.flag_2_probgood = sum(qc_matrix(:) == 2);
    counts.flag_3_suspect = sum(qc_matrix(:) == 3);
    counts.flag_4_bad = sum(qc_matrix(:) == 4);
    counts.flag_9_missing = sum(qc_matrix(:) == 9);
    counts.total = numel(qc_matrix);
    counts.pct_good = 100 * counts.flag_1_good / counts.total;
end

function print_qc_summary(summary)
% Print formatted QC summary
    vars = {'seapres', 'temp', 'cond', 'sal', 'sigma'};
    var_names = {'Sea Pressure', 'Temperature', 'Conductivity', 'Salinity', 'Density'};

    fprintf('  %-15s %7s %7s %7s %7s %7s\n', 'Variable', 'Good', 'Suspect', 'Bad', 'Missing', '%Good');
    fprintf('  %s\n', repmat('-', 1, 55));

    for i = 1:length(vars)
        c = summary.flags_count.(vars{i});
        fprintf('  %-15s %7d %7d %7d %7d %6.1f%%\n', ...
            var_names{i}, c.flag_1_good, c.flag_3_suspect, ...
            c.flag_4_bad, c.flag_9_missing, c.pct_good);
    end
end
