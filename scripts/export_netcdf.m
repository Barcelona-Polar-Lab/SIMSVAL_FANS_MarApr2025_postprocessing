function export_netcdf(matfile, ncfile, varargin)
% EXPORT_CTD_NETCDF  Convert processed CTD .mat to CF-1.8 / ACDD-1.3 NetCDF-4
%
% USAGE:
%   export_CTD('CTD_profiles.mat', 'CTD_profiles.nc')
%
% INPUT:
%   matfile  - path to .mat file
%   ncfile   - output NetCDF-4 file path
%
% OPTIONAL (name-value pairs):
%   'creator_name'        - default 'Nina Hoareau'
%   'creator_email'       - default 'nhoareau@icm.csic.es'
%   'creator_institution' - default 'ICM-CSIC, Barcelona'
%   'project'             - default 'ARICE-PONANT 2025'
%   'references'          - default ''
%   'overwrite'           - default true
%
% OUTPUT:
%   NetCDF-4 file with CF-1.8 conventions, ACDD-1.3 global attributes,
%   featureType = "profile", DSG profileIncomplete template.
%
% REQUIRES:
%   TEOS-10 GSW toolbox (for depth verification)
%
% Nina Hoareau, ICM-CSIC Barcelona, Jan. 2025

%% ========================================================================
%  PARSE INPUTS
%  ========================================================================
p = inputParser;
addRequired(p, 'matfile', @ischar);
addRequired(p, 'ncfile', @ischar);
addParameter(p, 'creator_name', 'Nina Hoareau', @ischar);
addParameter(p, 'creator_email', 'nhoareau@icm.csic.es', @ischar);
addParameter(p, 'creator_institution', 'ICM-CSIC, Barcelona', @ischar);
addParameter(p, 'project', 'ARICE 2025', @ischar);
addParameter(p, 'references', '', @ischar);
addParameter(p, 'overwrite', true, @islogical);
addParameter(p, 'sss_file', '', @ischar);
parse(p, matfile, ncfile, varargin{:});
opts = p.Results;

if exist(ncfile, 'file') && opts.overwrite
    delete(ncfile);
    fprintf('  Existing file deleted: %s\n', ncfile);
end

%% ========================================================================
%  LOAD DATA
%  ========================================================================
fprintf('\nLoading %s ...\n', matfile);
D = load(matfile);

[m, np] = size(D.ctdtemp);  % m = number of obs, np = number of profiles
fprintf('  %d num of obs x %d profiles\n', m, np);

% --- Check dimension consistency ---
assert(length(D.ctddate) == np, 'Dimension mismatch: ctddate (%d) vs ctdtemp (%d cols)', length(D.ctddate), np);
assert(length(D.ctdid) == np, 'Dimension mismatch: ctdid (%d) vs ctdtemp (%d cols)', length(D.ctdid), np);
assert(length(D.ctdlat) == np, 'Dimension mismatch: ctdlat (%d) vs ctdtemp (%d cols)', length(D.ctdlat), np);
assert(size(D.ctdseapres, 1) == m, 'Dimension mismatch: ctdseapres (%d) vs ctdtemp (%d rows)', size(D.ctdseapres,1), m);
fprintf('  [OK] Dimensions consistent\n');

%% ========================================================================
%  PREPARE DATA
%  ========================================================================

% --- Conductivity: mS/cm -> S/m ---
cond_Sm = D.ctdcond / 10;

% --- Time: MATLAB datetime -> POSIX seconds since 1970-01-01 ---
epoch = datetime(1970, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
dt_utc = D.ctddate(:)';
if isdatetime(dt_utc)
    if isempty(dt_utc.TimeZone)
        dt_utc.TimeZone = 'UTC';
    end
    time_posix = seconds(dt_utc - epoch);
else
    % If ctddate is datenum (legacy)
    dt_utc = datetime(D.ctddate(:)', 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');
    time_posix = seconds(dt_utc - epoch);
end

% --- Station IDs: convert to char matrix ---
if isstring(D.ctdid)
    id_cell = cellstr(D.ctdid);
elseif iscell(D.ctdid)
    id_cell = D.ctdid;
elseif ischar(D.ctdid)
    id_cell = cellstr(D.ctdid);
else
    id_cell = arrayfun(@(x) num2str(x), D.ctdid, 'UniformOutput', false);
end
max_strlen = max(cellfun(@length, id_cell));
station_char = char(id_cell);  % np x max_strlen

% --- Depth: use pre-computed 1-D vector from concat_oceanCast.m ---
depth_vec = D.ctddepth;  % (m x 1), already abs(gsw_z_from_p(...))
mean_lat = D.mean_lat;   % latitude used for depth calculation

% --- Geospatial/temporal bounds ---
lat_min = min(D.ctdlat, [], 'omitnan');
lat_max = max(D.ctdlat, [], 'omitnan');
lon_min = min(D.ctdlon, [], 'omitnan');
lon_max = max(D.ctdlon, [], 'omitnan');
pres_min = min(D.ctdseapres(:), [], 'omitnan');
pres_max = max(D.ctdseapres(:), [], 'omitnan');
time_start = min(dt_utc);
time_end   = max(dt_utc);

% --- n_samples: non-NaN observations per profile ---
n_samples_vec = int32(sum(~isnan(D.ctdseapres), 1));

% --- max_sea_pressure per profile ---
max_seapres_vec = max(D.ctdseapres, [], 1, 'omitnan');

% --- sampling_frequency per profile (Hz) ---
% Station C (IDs 'C_*') uses 2 Hz PONANT sensor; all others use 16 Hz ICM-CSIC.
sampling_freq_vec = zeros(1, np);
for ip = 1:np
    if startsWith(strtrim(id_cell{ip}), 'C_')
        sampling_freq_vec(ip) = 2;
    else
        sampling_freq_vec(ip) = 16;
    end
end

% --- Campaign: build char matrix (SIMSVAL / FANS) ---
campaign_cell = D.campaign;
campaign_strlen = max(cellfun(@length, campaign_cell));
campaign_char = char(campaign_cell);  % np x campaign_strlen

% --- Sample time: datenum -> POSIX seconds (NaN-safe) ---
sample_time_posix = NaN(m, np);
valid_st = ~isnan(D.ctdsampltime);
if any(valid_st(:))
    dt_samp = datetime(D.ctdsampltime(valid_st), 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');
    sample_time_posix(valid_st) = seconds(dt_samp - epoch);
end

% --- SSS water samples (optional, from AutoSal CSV) ---
sss_vec = NaN(1, np);
if ~isempty(opts.sss_file) && exist(opts.sss_file, 'file')
    sss_tbl = readtable(opts.sss_file, 'TextType', 'string');
    sss_map = containers.Map(string(sss_tbl.station_id), sss_tbl.sss_autosal);
    for ip = 1:np
        sid = strtrim(id_cell{ip});
        if isKey(sss_map, sid)
            sss_vec(ip) = sss_map(sid);
        end
    end
    fprintf('  SSS: %d values loaded from %s\n', sum(~isnan(sss_vec)), opts.sss_file);
elseif ~isempty(opts.sss_file)
    warning('export_netcdf:sssFileNotFound', 'SSS file not found: %s — sss filled with NaN', opts.sss_file);
end

%% ========================================================================
%  CREATE NETCDF FILE - DIMENSIONS & VARIABLES
%  ========================================================================
fprintf('\nCreating NetCDF-4 file: %s\n', ncfile);

FV     = NaN;        % FillValue for double
FV_qc  = int8(-1);   % FillValue for QC flags (outside 0-9 range)

% ---- 1. Coordinate: time(profile) ----
nccreate(ncfile, 'time', ...
    'Dimensions', {'profile', np}, ...
    'Datatype', 'double', ...
    'FillValue', FV, ...
    'Format', 'netcdf4');

% ---- 2. Coordinate: latitude(profile) ----
nccreate(ncfile, 'latitude', ...
    'Dimensions', {'profile', np}, ...
    'Datatype', 'double', ...
    'FillValue', FV);

% ---- 3. Coordinate: longitude(profile) ----
nccreate(ncfile, 'longitude', ...
    'Dimensions', {'profile', np}, ...
    'Datatype', 'double', ...
    'FillValue', FV);

% ---- 4. Coordinate: sea_pressure(z, profile) ----
nccreate(ncfile, 'sea_pressure', ...
    'Dimensions', {'obs', m, 'profile', np}, ...
    'Datatype', 'double', ...
    'FillValue', FV, ...
    'DeflateLevel', 4, ...
    'Shuffle', true);

% ---- 5. Coordinate: depth(z, profile) ----
nccreate(ncfile, 'depth', ...
    'Dimensions', {'obs', m, 'profile', np}, ...
    'Datatype', 'double', ...
    'FillValue', FV, ...
    'DeflateLevel', 4, ...
    'Shuffle', true);

% ---- 6. station_id(name_strlen, profile) ----
nccreate(ncfile, 'station_id', ...
    'Dimensions', {'name_strlen', max_strlen, 'profile', np}, ...
    'Datatype', 'char');

% ---- 7-11. Data variables (z, profile) ----
data_vars = {'temperature', 'salinity', 'conductivity', ...
             'sigma_theta', 'profile_velocity'};
for iv = 1:length(data_vars)
    nccreate(ncfile, data_vars{iv}, ...
        'Dimensions', {'obs', m, 'profile', np}, ...
        'Datatype', 'double', ...
        'FillValue', FV, ...
        'DeflateLevel', 4, ...
        'Shuffle', true);
end

% ---- 12-16. QC flag variables (z, profile), int8 ----
qc_vars = {'temperature_qc', 'salinity_qc', 'conductivity_qc', ...
            'sigma_theta_qc', 'sea_pressure_qc'};
for iv = 1:length(qc_vars)
    nccreate(ncfile, qc_vars{iv}, ...
        'Dimensions', {'obs', m, 'profile', np}, ...
        'Datatype', 'int8', ...
        'FillValue', FV_qc, ...
        'DeflateLevel', 4, ...
        'Shuffle', true);
end

% ---- 17-19. Soak auxiliary variables (profile) ----
soak_vars = {'soak_duration', 'soak_depth', 'soak_n_filtered'};
for iv = 1:length(soak_vars)
    nccreate(ncfile, soak_vars{iv}, ...
        'Dimensions', {'profile', np}, ...
        'Datatype', 'double', ...
        'FillValue', FV);
end

% ---- 20-24. Atmospheric variables (profile) ----
atm_vars = {'atm_pressure_1', 'atm_pressure_2', ...
             'wind_speed', 'wind_direction', 'relative_humidity'};
for iv = 1:length(atm_vars)
    nccreate(ncfile, atm_vars{iv}, ...
        'Dimensions', {'profile', np}, ...
        'Datatype', 'double', ...
        'FillValue', FV);
end

% ---- 25. n_samples(profile) ----
nccreate(ncfile, 'n_samples', ...
    'Dimensions', {'profile', np}, ...
    'Datatype', 'int32', ...
    'FillValue', int32(-1));

% ---- 26. max_sea_pressure(profile) ----
nccreate(ncfile, 'max_sea_pressure', ...
    'Dimensions', {'profile', np}, ...
    'Datatype', 'double', ...
    'FillValue', FV);

% ---- 27. sampling_frequency(profile) ----
nccreate(ncfile, 'sampling_frequency', ...
    'Dimensions', {'profile', np}, ...
    'Datatype', 'double', ...
    'FillValue', FV);

% ---- 28. sss(profile) ----
nccreate(ncfile, 'sss', ...
    'Dimensions', {'profile', np}, ...
    'Datatype', 'double', ...
    'FillValue', FV);

% ---- 29. campaign(profile, campaign_strlen) ----
nccreate(ncfile, 'campaign', ...
    'Dimensions', {'campaign_strlen', campaign_strlen, 'profile', np}, ...
    'Datatype', 'char');

% ---- 30. sample_time(obs, profile) ----
nccreate(ncfile, 'sample_time', ...
    'Dimensions', {'obs', m, 'profile', np}, ...
    'Datatype', 'double', ...
    'FillValue', FV, ...
    'DeflateLevel', 4, ...
    'Shuffle', true);

%% ========================================================================
%  WRITE VARIABLE ATTRIBUTES
%  ========================================================================
fprintf('  Writing variable attributes...\n');

% ============== COORDINATE VARIABLES ==============

% time
ncwriteatt(ncfile, 'time', 'standard_name', 'time');
ncwriteatt(ncfile, 'time', 'long_name', 'time of CTD profile');
ncwriteatt(ncfile, 'time', 'units', 'seconds since 1970-01-01T00:00:00Z');
ncwriteatt(ncfile, 'time', 'calendar', 'standard');
ncwriteatt(ncfile, 'time', 'axis', 'T');

% latitude
ncwriteatt(ncfile, 'latitude', 'standard_name', 'latitude');
ncwriteatt(ncfile, 'latitude', 'long_name', 'latitude of CTD profile');
ncwriteatt(ncfile, 'latitude', 'units', 'degrees_north');
ncwriteatt(ncfile, 'latitude', 'axis', 'Y');
ncwriteatt(ncfile, 'latitude', 'valid_min', -90.0);
ncwriteatt(ncfile, 'latitude', 'valid_max', 90.0);

% longitude
ncwriteatt(ncfile, 'longitude', 'standard_name', 'longitude');
ncwriteatt(ncfile, 'longitude', 'long_name', 'longitude of CTD profile');
ncwriteatt(ncfile, 'longitude', 'units', 'degrees_east');
ncwriteatt(ncfile, 'longitude', 'axis', 'X');
ncwriteatt(ncfile, 'longitude', 'valid_min', -180.0);
ncwriteatt(ncfile, 'longitude', 'valid_max', 180.0);

% sea_pressure
ncwriteatt(ncfile, 'sea_pressure', 'standard_name', 'sea_water_pressure');
ncwriteatt(ncfile, 'sea_pressure', 'long_name', 'sea water pressure (atmospheric pressure removed)');
ncwriteatt(ncfile, 'sea_pressure', 'units', 'dbar');
ncwriteatt(ncfile, 'sea_pressure', 'axis', 'Z');
ncwriteatt(ncfile, 'sea_pressure', 'positive', 'down');
ncwriteatt(ncfile, 'sea_pressure', 'ancillary_variables', 'sea_pressure_qc');

% depth
ncwriteatt(ncfile, 'depth', 'standard_name', 'depth');
ncwriteatt(ncfile, 'depth', 'long_name', 'depth below sea surface');
ncwriteatt(ncfile, 'depth', 'units', 'm');
ncwriteatt(ncfile, 'depth', 'positive', 'down');
ncwriteatt(ncfile, 'depth', 'comment', ...
    sprintf('Computed from sea_pressure using TEOS-10 gsw_z_from_p at mean latitude %.4f N', mean_lat));

% station_id
ncwriteatt(ncfile, 'station_id', 'long_name', 'station and profile identifier');
ncwriteatt(ncfile, 'station_id', 'cf_role', 'profile_id');

% ============== DATA VARIABLES ==============

% temperature
ncwriteatt(ncfile, 'temperature', 'standard_name', 'sea_water_temperature');
ncwriteatt(ncfile, 'temperature', 'long_name', 'sea water temperature from CTD');
ncwriteatt(ncfile, 'temperature', 'units', 'degree_Celsius');
ncwriteatt(ncfile, 'temperature', 'coordinates', 'time latitude longitude sea_pressure');
ncwriteatt(ncfile, 'temperature', 'ancillary_variables', 'temperature_qc');
ncwriteatt(ncfile, 'temperature', 'instrument', 'RBR Concerto CTD');

% salinity
ncwriteatt(ncfile, 'salinity', 'standard_name', 'sea_water_practical_salinity');
ncwriteatt(ncfile, 'salinity', 'long_name', 'practical salinity');
ncwriteatt(ncfile, 'salinity', 'units', '1');
ncwriteatt(ncfile, 'salinity', 'coordinates', 'time latitude longitude sea_pressure');
ncwriteatt(ncfile, 'salinity', 'ancillary_variables', 'salinity_qc');

% conductivity
ncwriteatt(ncfile, 'conductivity', 'standard_name', 'sea_water_electrical_conductivity');
ncwriteatt(ncfile, 'conductivity', 'long_name', 'sea water electrical conductivity');
ncwriteatt(ncfile, 'conductivity', 'units', 'S m-1');
ncwriteatt(ncfile, 'conductivity', 'coordinates', 'time latitude longitude sea_pressure');
ncwriteatt(ncfile, 'conductivity', 'ancillary_variables', 'conductivity_qc');
ncwriteatt(ncfile, 'conductivity', 'comment', 'Converted from mS/cm to S/m (divided by 10)');

% sigma_theta
ncwriteatt(ncfile, 'sigma_theta', 'standard_name', 'sea_water_sigma_theta');
ncwriteatt(ncfile, 'sigma_theta', 'long_name', 'potential density anomaly (sigma-theta)');
ncwriteatt(ncfile, 'sigma_theta', 'units', 'kg m-3');
ncwriteatt(ncfile, 'sigma_theta', 'coordinates', 'time latitude longitude sea_pressure');
ncwriteatt(ncfile, 'sigma_theta', 'ancillary_variables', 'sigma_theta_qc');

% profile_velocity
ncwriteatt(ncfile, 'profile_velocity', 'long_name', 'vertical profiling velocity during CTD downcast');
ncwriteatt(ncfile, 'profile_velocity', 'units', 'm s-1');
ncwriteatt(ncfile, 'profile_velocity', 'coordinates', 'time latitude longitude sea_pressure');
ncwriteatt(ncfile, 'profile_velocity', 'comment', 'Descent rate derived by RSKtools RSKderivevelocity');

% ============== QC FLAG VARIABLES ==============

qc_flag_values   = int8([0 1 2 3 4 9]);
qc_flag_meanings = 'no_qc_performed good_data probably_good_data suspect_or_doubtful_data bad_data missing_value';
qc_conventions   = 'SeaDataNet measurand qualifier flags';

qc_var_names  = {'temperature_qc', 'salinity_qc', 'conductivity_qc', ...
                 'sigma_theta_qc', 'sea_pressure_qc'};
qc_long_names = {'quality flag for sea water temperature', ...
                 'quality flag for practical salinity', ...
                 'quality flag for electrical conductivity', ...
                 'quality flag for sigma-theta', ...
                 'quality flag for sea water pressure'};

for iv = 1:length(qc_var_names)
    ncwriteatt(ncfile, qc_var_names{iv}, 'long_name', qc_long_names{iv});
    ncwriteatt(ncfile, qc_var_names{iv}, 'flag_values', qc_flag_values);
    ncwriteatt(ncfile, qc_var_names{iv}, 'flag_meanings', qc_flag_meanings);
    ncwriteatt(ncfile, qc_var_names{iv}, 'conventions', qc_conventions);
    ncwriteatt(ncfile, qc_var_names{iv}, 'valid_range', int8([0 9]));
    ncwriteatt(ncfile, qc_var_names{iv}, 'comment', ...
        ['Automated QC (RTQC-like, QARTOD-style). Tests applied: NaN check, gross range, ', ...
         'Arctic regional range, flat line, vertical gradient, density inversion, pressure monotonicity. ', ...
         'No delayed-mode validation performed.']);
end

% ============== SOAK AUXILIARY VARIABLES ==============

ncwriteatt(ncfile, 'soak_duration', 'long_name', 'duration of soak phase at start of profile');
ncwriteatt(ncfile, 'soak_duration', 'units', 's');
ncwriteatt(ncfile, 'soak_duration', 'comment', ...
    ['Soak removed using fixed_time method (default 20 s from profile start). ', ...
     'Exceptions: profiles 1A_1 and 7H_4 use 5 s (fast initial descent). ', ...
     'Function: RSKtrim_soak (v5).']);

ncwriteatt(ncfile, 'soak_depth', 'long_name', 'sea pressure at end of soak phase');
ncwriteatt(ncfile, 'soak_depth', 'units', 'dbar');
ncwriteatt(ncfile, 'soak_depth', 'comment', 'Sea pressure level where soak ends and valid descent begins');

ncwriteatt(ncfile, 'soak_n_filtered', 'long_name', 'number of data points removed during soak filtering');
ncwriteatt(ncfile, 'soak_n_filtered', 'units', '1');

% ============== ATMOSPHERIC VARIABLES ==============

ncwriteatt(ncfile, 'atm_pressure_1', 'long_name', 'atmospheric pressure from FerryBox sensor 1');
ncwriteatt(ncfile, 'atm_pressure_1', 'standard_name', 'air_pressure');
ncwriteatt(ncfile, 'atm_pressure_1', 'units', 'hPa');

ncwriteatt(ncfile, 'atm_pressure_2', 'long_name', 'atmospheric pressure from FerryBox sensor 2');
ncwriteatt(ncfile, 'atm_pressure_2', 'standard_name', 'air_pressure');
ncwriteatt(ncfile, 'atm_pressure_2', 'units', 'hPa');

ncwriteatt(ncfile, 'wind_speed', 'long_name', 'wind speed at time of profile');
ncwriteatt(ncfile, 'wind_speed', 'standard_name', 'wind_speed');
ncwriteatt(ncfile, 'wind_speed', 'units', 'm s-1');

ncwriteatt(ncfile, 'wind_direction', 'long_name', 'wind direction at time of profile');
ncwriteatt(ncfile, 'wind_direction', 'standard_name', 'wind_from_direction');
ncwriteatt(ncfile, 'wind_direction', 'units', 'degree');

ncwriteatt(ncfile, 'relative_humidity', 'long_name', 'relative humidity at time of profile');
ncwriteatt(ncfile, 'relative_humidity', 'standard_name', 'relative_humidity');
ncwriteatt(ncfile, 'relative_humidity', 'units', '%');

% ============== NEW PER-PROFILE VARIABLES ==============

ncwriteatt(ncfile, 'n_samples', 'long_name', 'number of valid observations per profile');
ncwriteatt(ncfile, 'n_samples', 'units', '1');
ncwriteatt(ncfile, 'n_samples', 'comment', 'Count of non-NaN sea_pressure values; excludes NaN-padding');

ncwriteatt(ncfile, 'max_sea_pressure', 'long_name', 'maximum sea water pressure reached during profile');
ncwriteatt(ncfile, 'max_sea_pressure', 'units', 'dbar');
ncwriteatt(ncfile, 'max_sea_pressure', 'comment', 'Maximum valid sea_pressure value per profile (NaN-padding excluded)');

ncwriteatt(ncfile, 'sampling_frequency', 'long_name', 'CTD sampling frequency');
ncwriteatt(ncfile, 'sampling_frequency', 'units', 'Hz');
ncwriteatt(ncfile, 'sampling_frequency', 'comment', ...
    ['16 Hz for ICM-CSIC RBR Concerto3 (S/N 237957); ', ...
     '2 Hz for PONANT RBR Concerto (station C). ', ...
     'Determined from station_id (C_* -> 2 Hz, all others -> 16 Hz).']);

ncwriteatt(ncfile, 'sss', 'long_name', 'sea surface salinity from AutoSal water sample');
ncwriteatt(ncfile, 'sss', 'standard_name', 'sea_water_practical_salinity');
ncwriteatt(ncfile, 'sss', 'units', '1');
ncwriteatt(ncfile, 'sss', 'comment', ...
    ['Surface salinity (0.2 m depth) from Guildline Autosal 8400B water samples ', ...
     'collected by zodiac. Accuracy: +/-0.001 PSU. ', ...
     'Complementary reference — not used for QC or calibration. ', ...
     'NaN where no water sample was collected.']);

% campaign
ncwriteatt(ncfile, 'campaign', 'long_name', 'campaign name');
ncwriteatt(ncfile, 'campaign', 'comment', 'SIMSVAL (March 2025) or FANS (April 2025)');

% sample_time
ncwriteatt(ncfile, 'sample_time', 'standard_name', 'time');
ncwriteatt(ncfile, 'sample_time', 'long_name', 'timestamp of each CTD sample');
ncwriteatt(ncfile, 'sample_time', 'units', 'seconds since 1970-01-01T00:00:00Z');
ncwriteatt(ncfile, 'sample_time', 'calendar', 'standard');
ncwriteatt(ncfile, 'sample_time', 'comment', ...
    'Per-sample UTC timestamp at native acquisition rate (16 Hz ICM-CSIC or 2 Hz PONANT). NaN for padding.');

%% ========================================================================
%  WRITE GLOBAL ATTRIBUTES (ACDD-1.3)
%  ========================================================================
fprintf('  Writing ACDD global attributes...\n');

% --- CF / ACDD ---
ncwriteatt(ncfile, '/', 'Conventions', 'CF-1.8, ACDD-1.3');
ncwriteatt(ncfile, '/', 'featureType', 'profile');
ncwriteatt(ncfile, '/', 'cdm_data_type', 'Profile');

% --- Title & Summary ---
ncwriteatt(ncfile, '/', 'title', ...
    'CTD profiles from ARICE-PONANT 2025 Arctic campaigns (SIMSVAL and FANS)');
ncwriteatt(ncfile, '/', 'summary', ...
    [sprintf('Processed CTD data from %d vertical profiles (SIMSVAL + FANS) ', np), ...
     'collected in Greenland fjords, Baffin Bay and adjacent Arctic waters ', ...
     'during the ARICE-PONANT 2025 expeditions aboard R/V Le Commandant Charcot. ', ...
     'Data were acquired with RBR Concerto3 (16 Hz, ICM-CSIC) and ', ...
     'RBR Concerto (2 Hz, PONANT) CTDs. Processing includes atmospheric ', ...
     'pressure correction, automatic soak detection and filtering, A2D correction, ', ...
     'despiking, CT lag adjustment (16 Hz only), smoothing, loop editing, ', ...
     'salinity and density derivation (TEOS-10). ', ...
     'Automated quality control (7 QARTOD/RTQC-style tests) was applied; ', ...
     'no delayed-mode validation (DMQC) was performed.']);

% --- Source & Platform ---
ncwriteatt(ncfile, '/', 'institution', opts.creator_institution);
ncwriteatt(ncfile, '/', 'source', ...
    'RBR Concerto3 CTD (16 Hz, S/N 237957, ICM-CSIC) and RBR Concerto CTD (2 Hz, PONANT)');
ncwriteatt(ncfile, '/', 'platform', 'R/V Le Commandant Charcot');
ncwriteatt(ncfile, '/', 'platform_vocabulary', 'https://vocab.nerc.ac.uk/collection/L06/current/');
ncwriteatt(ncfile, '/', 'instrument', 'RBR Concerto3 CTD, RBR Concerto CTD');
ncwriteatt(ncfile, '/', 'instrument_vocabulary', 'https://vocab.nerc.ac.uk/collection/L22/current/');

% --- Processing ---
ncwriteatt(ncfile, '/', 'processing_level', ...
    'L2: Post-processed with automated QC (RTQC-like, no delayed-mode validation)');
ncwriteatt(ncfile, '/', 'processing_software', 'RSKtools v3.6 (RBR Ltd.), TEOS-10 GSW, custom MATLAB scripts');
ncwriteatt(ncfile, '/', 'comment', ...
    ['Data at native temporal resolution (no bin averaging). Conductivity converted from mS/cm to S/m. ', ...
     'QC flags follow SeaDataNet convention (1=good, 2=probably good, 3=suspect, 4=bad, 9=missing). ', ...
     'Quality control consists of 7 automated tests (QARTOD/Argo RTQC-style): ', ...
     'NaN check, gross range, regional range, flat line, gradient, density inversion, pressure monotonicity. ', ...
     'No delayed-mode QC (DMQC) was performed: no comparison with climatologies or neighboring profiles, ', ...
     'no human expert validation, no sensor drift adjustment. ', ...
     'Soak phase removed using fixed_time method (20 s default; 5 s for 1A_1 and 7H_4). ', ...
     'Depth computed from sea_pressure at mean latitude using TEOS-10 gsw_z_from_p.']);

% --- Quality Control ---
ncwriteatt(ncfile, '/', 'qc_type', 'Automated (RTQC-like)');
ncwriteatt(ncfile, '/', 'qc_reference', 'QARTOD (IOOS), Argo QC Manual, SeaDataNet QC procedures');
ncwriteatt(ncfile, '/', 'qc_tests', ...
    'NaN check, gross range, regional range (Arctic), flat line, vertical gradient, density inversion, pressure monotonicity');
ncwriteatt(ncfile, '/', 'qc_note', ...
    ['Quality flags are based on automated tests only (Real-Time QC style). ', ...
     'No Delayed-Mode QC (DMQC) was performed: data were not compared against ', ...
     'climatologies, nearby Argo floats, or reference CTD profiles, and no ', ...
     'human expert validation was applied. Flag 1 indicates data passed all ', ...
     'automated tests but is not equivalent to DMQC-validated data.']);

% --- Creator ---
ncwriteatt(ncfile, '/', 'creator_name', opts.creator_name);
ncwriteatt(ncfile, '/', 'creator_email', opts.creator_email);
ncwriteatt(ncfile, '/', 'creator_institution', opts.creator_institution);
ncwriteatt(ncfile, '/', 'creator_role', 'Data processing and quality control');

% --- Principal Investigators ---
ncwriteatt(ncfile, '/', 'pi_name', 'Carolina Gabarro, Marta Umbert');
ncwriteatt(ncfile, '/', 'pi_email', 'cgabarro@icm.csic.es, mumbert@icm.csic.es');
ncwriteatt(ncfile, '/', 'pi_institution', 'ICM-CSIC, Barcelona');

% --- Contributors ---
ncwriteatt(ncfile, '/', 'contributor_name', 'Nina Hoareau, Carolina Gabarro, Marta Umbert');
ncwriteatt(ncfile, '/', 'contributor_role', 'Data processor, Principal Investigator, Principal Investigator');

% --- Project ---
ncwriteatt(ncfile, '/', 'project', opts.project);
if ~isempty(opts.references)
    ncwriteatt(ncfile, '/', 'references', opts.references);
end

% --- Geospatial Bounds ---
ncwriteatt(ncfile, '/', 'geospatial_lat_min', lat_min);
ncwriteatt(ncfile, '/', 'geospatial_lat_max', lat_max);
ncwriteatt(ncfile, '/', 'geospatial_lat_units', 'degrees_north');
ncwriteatt(ncfile, '/', 'geospatial_lon_min', lon_min);
ncwriteatt(ncfile, '/', 'geospatial_lon_max', lon_max);
ncwriteatt(ncfile, '/', 'geospatial_lon_units', 'degrees_east');
ncwriteatt(ncfile, '/', 'geospatial_vertical_min', pres_min);
ncwriteatt(ncfile, '/', 'geospatial_vertical_max', pres_max);
ncwriteatt(ncfile, '/', 'geospatial_vertical_units', 'dbar');
ncwriteatt(ncfile, '/', 'geospatial_vertical_positive', 'down');

% --- Time Coverage ---
ncwriteatt(ncfile, '/', 'time_coverage_start', datestr(time_start, 'yyyy-mm-ddTHH:MM:SSZ'));
ncwriteatt(ncfile, '/', 'time_coverage_end',   datestr(time_end,   'yyyy-mm-ddTHH:MM:SSZ'));

% --- History & Dates ---
ncwriteatt(ncfile, '/', 'history', ...
    [datestr(now, 'yyyy-mm-ddTHH:MM:SSZ'), ' Created by export_CTD_netcdf.m from ', matfile]);
ncwriteatt(ncfile, '/', 'date_created',  datestr(now, 'yyyy-mm-ddTHH:MM:SSZ'));
ncwriteatt(ncfile, '/', 'date_modified', datestr(now, 'yyyy-mm-ddTHH:MM:SSZ'));

%% ========================================================================
%  WRITE DATA
%  ========================================================================
fprintf('  Writing data...\n');

% --- Coordinates ---
ncwrite(ncfile, 'time', time_posix);
ncwrite(ncfile, 'latitude', D.ctdlat(:)');
ncwrite(ncfile, 'longitude', D.ctdlon(:)');
ncwrite(ncfile, 'sea_pressure', D.ctdseapres);
ncwrite(ncfile, 'depth', depth_vec);

% --- Station IDs (char array, transposed for NetCDF) ---
ncwrite(ncfile, 'station_id', station_char');

% --- 2-D data variables ---
ncwrite(ncfile, 'temperature',      D.ctdtemp);
ncwrite(ncfile, 'salinity',         D.ctdsal);
ncwrite(ncfile, 'conductivity',     cond_Sm);
ncwrite(ncfile, 'sigma_theta',      D.ctdsigma);
ncwrite(ncfile, 'profile_velocity', D.ctdvelprof);

% --- QC flags (int8) ---
ncwrite(ncfile, 'temperature_qc',   D.ctdtemp_qc);
ncwrite(ncfile, 'salinity_qc',      D.ctdsal_qc);
ncwrite(ncfile, 'conductivity_qc',  D.ctdcond_qc);
ncwrite(ncfile, 'sigma_theta_qc',   D.ctdsigma_qc);
ncwrite(ncfile, 'sea_pressure_qc',  D.ctdseapres_qc);

% --- Soak info ---
ncwrite(ncfile, 'soak_duration',    D.soak_duration_s(:)');
ncwrite(ncfile, 'soak_depth',       D.soak_depth_dbar(:)');
ncwrite(ncfile, 'soak_n_filtered',  double(D.soak_n_filtered(:)'));

% --- Atmospheric ---
ncwrite(ncfile, 'atm_pressure_1',   D.atmhPa1(:)');
ncwrite(ncfile, 'atm_pressure_2',   D.atmhPa2(:)');
ncwrite(ncfile, 'wind_speed',       D.atmwindspeed(:)');
ncwrite(ncfile, 'wind_direction',   D.atmwinddirection(:)');
ncwrite(ncfile, 'relative_humidity',D.atmrelativehumidity(:)');

% --- New per-profile variables ---
ncwrite(ncfile, 'n_samples',          n_samples_vec);
ncwrite(ncfile, 'max_sea_pressure',   max_seapres_vec);
ncwrite(ncfile, 'sampling_frequency', sampling_freq_vec);
ncwrite(ncfile, 'sss',                sss_vec);
ncwrite(ncfile, 'campaign',           campaign_char');
ncwrite(ncfile, 'sample_time',        sample_time_posix);

%% ========================================================================
%  VERIFICATION
%  ========================================================================
fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('  EXPORT NETCDF - VERIFICATION\n');
fprintf('════════════════════════════════════════════════════════════════\n');

info = ncinfo(ncfile);
fprintf('  File:   %s\n', ncfile);
fprintf('  Format: %s\n', info.Format);

fprintf('\n  Dimensions:\n');
for id = 1:length(info.Dimensions)
    fprintf('    %-12s = %d\n', info.Dimensions(id).Name, info.Dimensions(id).Length);
end

fprintf('\n  Variables (%d):\n', length(info.Variables));
for iv = 1:length(info.Variables)
    v = info.Variables(iv);
    dims_str = strjoin({v.Dimensions.Name}, ' x ');
    fprintf('    %-22s  (%s)\n', v.Name, dims_str);
end

% Roundtrip check
temp_check = ncread(ncfile, 'temperature');
max_diff_t = max(abs(temp_check(:) - D.ctdtemp(:)), [], 'omitnan');
fprintf('\n  Roundtrip temperature: max diff = %.2e\n', max_diff_t);

cond_check = ncread(ncfile, 'conductivity');
max_diff_c = max(abs(cond_check(:) - cond_Sm(:)), [], 'omitnan');
fprintf('  Roundtrip conductivity (S/m): max diff = %.2e\n', max_diff_c);

tqc = ncread(ncfile, 'temperature_qc');
fprintf('  Temperature QC - unique values: %s\n', mat2str(unique(tqc(:))'));

sid = ncread(ncfile, 'station_id');
fprintf('  First station IDs: ');
for k = 1:min(5, np)
    fprintf('%s  ', strtrim(sid(:,k)'));
end
fprintf('\n');

t_check = ncread(ncfile, 'time');
fprintf('  Period: %s to %s\n', ...
    datestr(datetime(min(t_check), 'ConvertFrom', 'posixtime')), ...
    datestr(datetime(max(t_check), 'ConvertFrom', 'posixtime')));

d = dir(ncfile);
fprintf('  File size: %.1f KB\n', d.bytes/1024);

fprintf('\n  For full inspection:\n');
fprintf('    ncdisp(''%s'')\n', ncfile);
fprintf('  For CF validation:\n');
fprintf('    cfchecks %s\n', ncfile);
fprintf('════════════════════════════════════════════════════════════════\n');

end
