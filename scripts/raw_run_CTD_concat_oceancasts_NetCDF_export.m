%% raw_run_CTD_concat_oceancasts_NetCDF_export.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nina HOAREAU, Avril 2026 (BPL, ICM-CSIC)
% ---
% Script for exporting raw CTD data to NetCDF-4 (CF-1.8 / ACDD-1.3).
% version 5
%
% PREREQUISITES:
%   1. Run raw_concat_SIMSVAL_matfile.m
%         → RAW_CTD_SIMSVAL_oceanCasts.mat
%   2. Run raw_concat_FANS_matfile.m
%         → RAW_CTD_FANS_oceanCasts.mat
%
% Structure of the output NetCDF file:
%   Dimensions:
%     profile (80)          — one per CTD profile
%     obs     (obs_max)     — max samples across all profiles (NaN-padded)
%     name_strlen           — max length of station IDs
%
%   Per-sample variables (obs × profile):
%     sample_time           — absolute timestamp (POSIX s since 1970-01-01)
%     pressure              — raw absolute pressure (dbar)
%     sea_pressure          — sea pressure = pressure - Patm (dbar)
%     temperature           — raw temperature (°C)
%     conductivity          — raw conductivity (mS/cm)
%
%   Per-profile variables (profile):
%     time                  — profile start time (POSIX s, CF coordinate)
%     latitude, longitude   — geographic coordinates
%     station_id            — station identifier (cf_role = profile_id)
%     n_samples             — number of valid samples (rest = NaN padding)
%     sampling_frequency    — sensor frequency (Hz: 16 or 2)
%     max_pressure          — maximum pressure reached in profile (dbar)
%     sss                   — near-surface salinity from water sample (PSU, complementary var.)
%                             Source: ancillary_data/SSS_autosal.csv
%     atm_pressure_1/2      — FerryBox atmospheric pressure (hPa)
%     wind_speed            — wind speed (m/s)
%     wind_direction        — wind direction (degrees)
%     relative_humidity     — relative humidity (%)
%
% Output: RAW_ARICE-2025_Greenland_CTD.nc
%
% NOTE ON PROCESSING LEVEL:
%   These data are L0: no processing has been applied. The data contain
%   the surface soak phase, pressure loops, potential spikes, and the
%   atmospheric pressure offset.
%   Use the PROCESSED file for scientific analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% ════════════════════════════════════════════════════════════════════════
%  PATHS
%  ════════════════════════════════════════════════════════════════════════

out_path = '../outputs/';
matfile_S = [out_path 'RAW_CTD_SIMSVAL_oceanCasts.mat'];
matfile_F = [out_path 'RAW_CTD_FANS_oceanCasts.mat'];
ncfile    = [out_path 'RAW_ARICE-2025_Greenland_CTD.nc'];

disp(' ')
disp('╔════════════════════════════════════════════════════════════════════╗')
disp('║  EXPORT NETCDF RAW — CTD ARICE 2025 Greenland                      ║')
disp('╚════════════════════════════════════════════════════════════════════╝')

%% ════════════════════════════════════════════════════════════════════════
%  STEP 1: LOAD .MAT FILES
%  ════════════════════════════════════════════════════════════════════════
disp(' ')
disp('--- Step 1: Load .mat files ---')
fprintf('  SIMSVAL: %s\n', matfile_S)
S = load(matfile_S);
fprintf('  FANS:    %s\n', matfile_F)
F = load(matfile_F);

np_S = length(S.ctdid);
np_F = length(F.ctdid);
np   = np_S + np_F;
fprintf('  %d profiles SIMSVAL + %d profiles FANS = %d profiles total\n', np_S, np_F, np)

%% ════════════════════════════════════════════════════════════════════════
%  STEP 2: CONCATENATION
%  ════════════════════════════════════════════════════════════════════════
disp(' ')
disp('--- Step 2: Concatenation ---')

raw_T   = [S.raw_T   F.raw_T];
raw_C   = [S.raw_C   F.raw_C];
raw_P   = [S.raw_P   F.raw_P];
raw_t   = [S.raw_t   F.raw_t];
n_samp  = [S.n_samp  F.n_samp];
freq    = [S.freq    F.freq];
ctdid   = [S.ctdid   F.ctdid];
ctdlat  = [S.ctdlat(:)'  F.ctdlat(:)'];
ctdlon  = [S.ctdlon(:)'  F.ctdlon(:)'];
ctddate = [S.ctddate(:)' F.ctddate(:)'];

atmhPa1          = [S.atmhPa1(:)'          F.atmhPa1(:)'];
atmhPa2          = [S.atmhPa2(:)'          F.atmhPa2(:)'];
atmwindspeed     = [S.atmwindspeed(:)'     F.atmwindspeed(:)'];
atmwinddirection = [S.atmwinddirection(:)' F.atmwinddirection(:)'];
atmrelativehumidity = [S.atmrelativehumidity(:)' F.atmrelativehumidity(:)'];

obs_max = max(n_samp);
fprintf('  obs_max (max samples per profile) = %d\n', obs_max)
fprintf('  Frequencies: %d Hz min, %d Hz max\n', min(freq), max(freq))

%% ════════════════════════════════════════════════════════════════════════
%  STEP 3: 2D PADDING (obs_max × np)
%  ════════════════════════════════════════════════════════════════════════
disp(' ')
disp('--- Step 3: 2D padding (obs_max × profile) ---')

T_2d  = NaN(obs_max, np);
C_2d  = NaN(obs_max, np);
P_2d  = NaN(obs_max, np);
SP_2d = NaN(obs_max, np);   % Sea pressure = P - Patm
t_2d  = NaN(obs_max, np);   % Sample timestamps (POSIX s)

for ip = 1:np
    ns = n_samp(ip);
    T_2d(1:ns, ip)  = raw_T{ip};
    C_2d(1:ns, ip)  = raw_C{ip};
    P_2d(1:ns, ip)  = raw_P{ip};
    t_2d(1:ns, ip)  = raw_t{ip};

    % Sea pressure: P - Patm [dbar]
    % Patm = mean of the 2 FerryBox sensors, converted hPa → dbar (÷100)
    Patm_ip = ((atmhPa1(ip) + atmhPa2(ip)) / 2) / 100;
    if isnan(Patm_ip)
        Patm_ip = 10.1325;   % standard atmosphere if FerryBox data absent
    end
    SP_2d(1:ns, ip) = raw_P{ip} - Patm_ip;
end

% Maximum pressure reached per profile
max_pres = NaN(1, np);
for ip = 1:np
    max_pres(ip) = max(raw_P{ip}, [], 'omitnan');
end

fprintf('  Final dimensions: %d × %d (obs × profiles)\n', obs_max, np)
fprintf('  Global max pressure: %.1f dbar\n', max(max_pres, [], 'omitnan'))

%% ════════════════════════════════════════════════════════════════════════
%  STEP 4: NEAR-SURFACE SALINITY FROM WATER SAMPLES (complementary variable)
%  Source: ancillary_data/SSS_autosal.csv
%  NaN = no water sample available for this profile
%  ════════════════════════════════════════════════════════════════════════
disp(' ')
disp('--- Step 4: Near-surface salinity AutoSal ---')

sss_file = '../ancillary_data/SSS_autosal.csv';
sss_t    = readtable(sss_file);
sss_map  = containers.Map(sss_t.station_id, sss_t.sss_autosal);

sss = NaN(1, np);
for ip = 1:np
    if isKey(sss_map, ctdid{ip})
        sss(ip) = sss_map(ctdid{ip});
    end
end
fprintf('  SSS: %d valid values / %d profiles\n', sum(~isnan(sss)), np)

%% ════════════════════════════════════════════════════════════════════════
%  STEP 5: PREPARE NETCDF VARIABLES
%  ════════════════════════════════════════════════════════════════════════
disp(' ')
disp('--- Step 5: Prepare NetCDF variables ---')

epoch = datetime(1970, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');

% Profile start time (POSIX s)
time_profile = NaN(1, np);
for ip = 1:np
    if ~isnat(ctddate(ip))
        time_profile(ip) = seconds(ctddate(ip) - epoch);
    end
end

% Station IDs → char matrix
max_strlen = max(cellfun(@length, ctdid));
station_char = char(ctdid);   % np × max_strlen

% Geospatial bounds
lat_min = min(ctdlat, [], 'omitnan');
lat_max = max(ctdlat, [], 'omitnan');
lon_min = min(ctdlon, [], 'omitnan');
lon_max = max(ctdlon, [], 'omitnan');
pres_min = min(SP_2d(:), [], 'omitnan');
pres_max = max(SP_2d(:), [], 'omitnan');
time_start = min(ctddate(~isnat(ctddate)));
time_end   = max(ctddate(~isnat(ctddate)));

fprintf('  obs_max = %d, np = %d, name_strlen = %d\n', obs_max, np, max_strlen)

%% ════════════════════════════════════════════════════════════════════════
%  STEP 6: CREATE NETCDF-4 FILE
%  ════════════════════════════════════════════════════════════════════════
disp(' ')
disp('--- Step 6: Create NetCDF-4 ---')

if exist(ncfile, 'file')
    delete(ncfile);
    fprintf('  Existing file deleted.\n')
end

FV = NaN;   % FillValue for double

% ---- DIMENSIONS AND VARIABLES ----

% 1. time(profile) — CF coordinate
nccreate(ncfile, 'time', ...
    'Dimensions', {'profile', np}, ...
    'Datatype', 'double', ...
    'FillValue', FV, ...
    'Format', 'netcdf4');

% 2. latitude(profile)
nccreate(ncfile, 'latitude', ...
    'Dimensions', {'profile', np}, ...
    'Datatype', 'double', ...
    'FillValue', FV);

% 3. longitude(profile)
nccreate(ncfile, 'longitude', ...
    'Dimensions', {'profile', np}, ...
    'Datatype', 'double', ...
    'FillValue', FV);

% 4. station_id(name_strlen, profile)
nccreate(ncfile, 'station_id', ...
    'Dimensions', {'name_strlen', max_strlen, 'profile', np}, ...
    'Datatype', 'char');

% 5. n_samples(profile)
nccreate(ncfile, 'n_samples', ...
    'Dimensions', {'profile', np}, ...
    'Datatype', 'int32', ...
    'FillValue', int32(-1));

% 6. sampling_frequency(profile)
nccreate(ncfile, 'sampling_frequency', ...
    'Dimensions', {'profile', np}, ...
    'Datatype', 'double', ...
    'FillValue', FV);

% 7. max_pressure(profile)
nccreate(ncfile, 'max_pressure', ...
    'Dimensions', {'profile', np}, ...
    'Datatype', 'double', ...
    'FillValue', FV);

% 8. sss(profile)
nccreate(ncfile, 'sss', ...
    'Dimensions', {'profile', np}, ...
    'Datatype', 'double', ...
    'FillValue', FV);

% 9-10. Atmospheric pressure (profile)
for vname = {'atm_pressure_1', 'atm_pressure_2', 'wind_speed', 'wind_direction', 'relative_humidity'}
    nccreate(ncfile, vname{1}, ...
        'Dimensions', {'profile', np}, ...
        'Datatype', 'double', ...
        'FillValue', FV);
end

% 11. sample_time(obs, profile) — timestamp of each raw sample
nccreate(ncfile, 'sample_time', ...
    'Dimensions', {'obs', obs_max, 'profile', np}, ...
    'Datatype', 'double', ...
    'FillValue', FV, ...
    'DeflateLevel', 4, ...
    'Shuffle', true);

% 12. pressure(obs, profile) — raw absolute pressure
nccreate(ncfile, 'pressure', ...
    'Dimensions', {'obs', obs_max, 'profile', np}, ...
    'Datatype', 'double', ...
    'FillValue', FV, ...
    'DeflateLevel', 4, ...
    'Shuffle', true);

% 13. sea_pressure(obs, profile) — raw sea pressure
nccreate(ncfile, 'sea_pressure', ...
    'Dimensions', {'obs', obs_max, 'profile', np}, ...
    'Datatype', 'double', ...
    'FillValue', FV, ...
    'DeflateLevel', 4, ...
    'Shuffle', true);

% 14. temperature(obs, profile)
nccreate(ncfile, 'temperature', ...
    'Dimensions', {'obs', obs_max, 'profile', np}, ...
    'Datatype', 'double', ...
    'FillValue', FV, ...
    'DeflateLevel', 4, ...
    'Shuffle', true);

% 15. conductivity(obs, profile)
nccreate(ncfile, 'conductivity', ...
    'Dimensions', {'obs', obs_max, 'profile', np}, ...
    'Datatype', 'double', ...
    'FillValue', FV, ...
    'DeflateLevel', 4, ...
    'Shuffle', true);

%% ════════════════════════════════════════════════════════════════════════
%  STEP 7: VARIABLE ATTRIBUTES
%  ════════════════════════════════════════════════════════════════════════
disp('  Writing variable attributes...')

% time
ncwriteatt(ncfile, 'time', 'standard_name', 'time');
ncwriteatt(ncfile, 'time', 'long_name', 'time of first sample in CTD profile (profile start)');
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
ncwriteatt(ncfile, 'latitude', 'comment', 'NaN for profiles with unknown position (7 profiles)');

% longitude
ncwriteatt(ncfile, 'longitude', 'standard_name', 'longitude');
ncwriteatt(ncfile, 'longitude', 'long_name', 'longitude of CTD profile');
ncwriteatt(ncfile, 'longitude', 'units', 'degrees_east');
ncwriteatt(ncfile, 'longitude', 'axis', 'X');
ncwriteatt(ncfile, 'longitude', 'valid_min', -180.0);
ncwriteatt(ncfile, 'longitude', 'valid_max', 180.0);

% station_id
ncwriteatt(ncfile, 'station_id', 'long_name', 'station and profile identifier');
ncwriteatt(ncfile, 'station_id', 'cf_role', 'profile_id');

% n_samples
ncwriteatt(ncfile, 'n_samples', 'long_name', 'number of valid raw samples in this profile');
ncwriteatt(ncfile, 'n_samples', 'units', '1');
ncwriteatt(ncfile, 'n_samples', 'comment', ...
    'Samples 1:n_samples are valid; samples n_samples+1:obs_max are NaN (padding).');

% sampling_frequency
ncwriteatt(ncfile, 'sampling_frequency', 'long_name', 'CTD sampling frequency');
ncwriteatt(ncfile, 'sampling_frequency', 'units', 'Hz');
ncwriteatt(ncfile, 'sampling_frequency', 'comment', ...
    '16 Hz RBR Concerto3 CTD; 2 Hz RBR Concerto3 CTD');

% max_pressure
ncwriteatt(ncfile, 'max_pressure', 'long_name', 'maximum absolute pressure reached during profile');
ncwriteatt(ncfile, 'max_pressure', 'units', 'dbar');

% sss
ncwriteatt(ncfile, 'sss', 'standard_name', 'sea_water_practical_salinity');
ncwriteatt(ncfile, 'sss', 'long_name', 'near-surface practical salinity from water samples analized with Autosal');
ncwriteatt(ncfile, 'sss', 'units', '1');
ncwriteatt(ncfile, 'sss', 'comment', ...
    ['Complementary variable: salinity measured from water sample, ', ...
     'collected at the CTD station. AutoSal analysis. ', ...
     'NaN = no Niskin sample available for this profile. ', ...
     'NOT used as a calibration reference (see dataset documentation for justification).']);

% sample_time
ncwriteatt(ncfile, 'sample_time', 'long_name', 'absolute timestamp of each raw CTD sample');
ncwriteatt(ncfile, 'sample_time', 'units', 'seconds since 1970-01-01T00:00:00Z');
ncwriteatt(ncfile, 'sample_time', 'calendar', 'standard');
ncwriteatt(ncfile, 'sample_time', 'comment', ...
    'Per-sample timestamps. Valid for indices 1:n_samples; NaN for padding beyond n_samples.');

% pressure
ncwriteatt(ncfile, 'pressure', 'long_name', 'raw absolute pressure from CTD sensor');
ncwriteatt(ncfile, 'pressure', 'standard_name', 'pressure');
ncwriteatt(ncfile, 'pressure', 'units', 'dbar');
ncwriteatt(ncfile, 'pressure', 'coordinates', 'sample_time latitude longitude');
ncwriteatt(ncfile, 'pressure', 'comment', ...
    ['Raw absolute pressure as recorded by the RBR pressure sensor. ', ...
     'Includes atmospheric pressure contribution (~10.1 dbar at surface). ', ...
     'To obtain sea water pressure, subtract Patm = (atm_pressure_1 + atm_pressure_2) / 2 / 100 dbar, ', ...
     'or use the pre-computed sea_pressure variable.']);

% sea_pressure
ncwriteatt(ncfile, 'sea_pressure', 'long_name', 'sea water pressure (atmospheric pressure removed)');
ncwriteatt(ncfile, 'sea_pressure', 'standard_name', 'sea_water_pressure');
ncwriteatt(ncfile, 'sea_pressure', 'units', 'dbar');
ncwriteatt(ncfile, 'sea_pressure', 'coordinates', 'sample_time latitude longitude');
ncwriteatt(ncfile, 'sea_pressure', 'positive', 'down');
ncwriteatt(ncfile, 'sea_pressure', 'comment', ...
    ['Sea pressure = pressure - Patm, where Patm = (atm_pressure_1 + atm_pressure_2) / 2 / 100 dbar. ', ...
     'Atmospheric pressure from FerryBox boat, hourly means matched by nearest neighbour to profile start time. ', ...
     'If FerryBox data unavailable (NaN), Patm = 10.1325 dbar (standard atmosphere) was used.']);

% temperature
ncwriteatt(ncfile, 'temperature', 'standard_name', 'sea_water_temperature');
ncwriteatt(ncfile, 'temperature', 'long_name', 'raw sea water temperature from CTD');
ncwriteatt(ncfile, 'temperature', 'units', 'degree_Celsius');
ncwriteatt(ncfile, 'temperature', 'coordinates', 'sample_time latitude longitude sea_pressure');
ncwriteatt(ncfile, 'temperature', 'instrument', 'RBR Concerto CTD');
ncwriteatt(ncfile, 'temperature', 'comment', ...
    ['Raw temperature. No processing applied (no soak filtering, no despiking, ', ...
     'no CT lag correction, no smoothing, no loop editing). ', ...
     'Contains soak phase at start of profile and potential spikes.']);

% conductivity
ncwriteatt(ncfile, 'conductivity', 'standard_name', 'sea_water_electrical_conductivity');
ncwriteatt(ncfile, 'conductivity', 'long_name', 'raw sea water electrical conductivity from CTD');
ncwriteatt(ncfile, 'conductivity', 'units', 'mS cm-1');
ncwriteatt(ncfile, 'conductivity', 'coordinates', 'sample_time latitude longitude sea_pressure');
ncwriteatt(ncfile, 'conductivity', 'comment', ...
    ['Raw conductivity'...
     'No processing applied. No CT lag correction, no hold correction.']);

% Atmospheric variables
ncwriteatt(ncfile, 'atm_pressure_1', 'long_name', 'atmospheric pressure from FerryBox sensor 1');
ncwriteatt(ncfile, 'atm_pressure_1', 'standard_name', 'air_pressure');
ncwriteatt(ncfile, 'atm_pressure_1', 'units', 'hPa');
ncwriteatt(ncfile, 'atm_pressure_1', 'comment', 'Hourly mean, matched to profile start time by nearest neighbour');

ncwriteatt(ncfile, 'atm_pressure_2', 'long_name', 'atmospheric pressure from FerryBox sensor 2');
ncwriteatt(ncfile, 'atm_pressure_2', 'standard_name', 'air_pressure');
ncwriteatt(ncfile, 'atm_pressure_2', 'units', 'hPa');
ncwriteatt(ncfile, 'atm_pressure_2', 'comment', 'Hourly mean, matched to profile start time by nearest neighbour');

ncwriteatt(ncfile, 'wind_speed', 'long_name', 'wind speed at time of profile');
ncwriteatt(ncfile, 'wind_speed', 'standard_name', 'wind_speed');
ncwriteatt(ncfile, 'wind_speed', 'units', 'm s-1');

ncwriteatt(ncfile, 'wind_direction', 'long_name', 'wind direction at time of profile');
ncwriteatt(ncfile, 'wind_direction', 'standard_name', 'wind_from_direction');
ncwriteatt(ncfile, 'wind_direction', 'units', 'degree');

ncwriteatt(ncfile, 'relative_humidity', 'long_name', 'relative humidity at time of profile');
ncwriteatt(ncfile, 'relative_humidity', 'standard_name', 'relative_humidity');
ncwriteatt(ncfile, 'relative_humidity', 'units', '%');

%% ════════════════════════════════════════════════════════════════════════
%  STEP 8: GLOBAL ATTRIBUTES (ACDD-1.3 / CF-1.8)
%  ════════════════════════════════════════════════════════════════════════
disp('  Writing global attributes...')

ncwriteatt(ncfile, '/', 'Conventions', 'CF-1.8, ACDD-1.3');
ncwriteatt(ncfile, '/', 'featureType', 'profile');
ncwriteatt(ncfile, '/', 'cdm_data_type', 'Profile');

ncwriteatt(ncfile, '/', 'title', ...
    'Raw CTD profiles from ARICE 2025 Greenland campaigns (SIMSVAL and FANS)');

ncwriteatt(ncfile, '/', 'summary', ...
    [sprintf('Raw (unprocessed) CTD data from %d vertical profiles collected during ', np), ...
     'the SIMSVAL and FANS campaigns of the ARICE 2025 Greenland expeditions aboard ', ...
     'R/V Le Commandant Charcot in western Greenland fjords. ', ...
     'Data were acquired with a RBR Concerto3 (16 Hz, S/N 237957) and ', ...
     'a RBR Concerto (2 Hz, S/N 237329) CTD. ', ...
     'This file contains L0 data: no soak filtering, no despiking, no CT lag correction, ', ...
     'no smoothing, no loop editing, no bin averaging. ', ...
     'The only metadata correction applied is the computation of sea_pressure = pressure - Patm. ', ...
     'For science-ready data, use the L2 processed file: ARICE_2025_Greenland_CTD.nc']);

ncwriteatt(ncfile, '/', 'institution', 'ICM-CSIC, Barcelona');
ncwriteatt(ncfile, '/', 'source', ...
    'RBR Concerto3 CTD (16 Hz, S/N 237957) and RBR Concerto CTD (2 Hz, S/N 237329)');
ncwriteatt(ncfile, '/', 'platform', 'R/V Le Commandant Charcot');
ncwriteatt(ncfile, '/', 'instrument', 'RBR Concerto3 CTD');
ncwriteatt(ncfile, '/', 'instrument_vocabulary', 'https://vocab.nerc.ac.uk/collection/L22/current/');

ncwriteatt(ncfile, '/', 'processing_level', ...
    'L0: Raw data, no processing applied. Sea pressure computed from raw pressure and FerryBox Atmospheric Pressure only.');
ncwriteatt(ncfile, '/', 'processing_software', 'RSKtools v3.6 (RBR Ltd.), custom MATLAB scripts');

ncwriteatt(ncfile, '/', 'comment', ...
    ['2D structure: dimensions (obs x profile). ', ...
     'obs = sample index (max samples across all profiles, NaN-padded). ', ...
     'n_samples(profile) indicates data length for each profile. ', ...
     'No quality control flags applied to this L0 file. ', ...
     'Profiles include surface soak phase and may contain spikes and pressure loops. ', ...
     '7 profiles have NaN coordinates (position unknown at time of cast). ', ...
     'See companion post-processed file ARICE_2025_Greenland_CTD.nc for processed data with QC flags.']);

ncwriteatt(ncfile, '/', 'creator_name', 'Nina Hoareau');
ncwriteatt(ncfile, '/', 'creator_email', 'nhoareau@icm.csic.es');
ncwriteatt(ncfile, '/', 'creator_institution', 'BPL (ICM-CSIC), Barcelona');
ncwriteatt(ncfile, '/', 'creator_role', 'Data processing');

ncwriteatt(ncfile, '/', 'pi_name', 'Marta Umbert, Carolina Gabarro');
ncwriteatt(ncfile, '/', 'pi_email', 'mumbert@icm.csic.es, cgabarro@icm.csic.es');
ncwriteatt(ncfile, '/', 'pi_institution', 'BPL (ICM-CSIC), Barcelona');

ncwriteatt(ncfile, '/', 'project', 'ARICE 2025');

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

ncwriteatt(ncfile, '/', 'time_coverage_start', datestr(time_start, 'yyyy-mm-ddTHH:MM:SSZ'));
ncwriteatt(ncfile, '/', 'time_coverage_end',   datestr(time_end,   'yyyy-mm-ddTHH:MM:SSZ'));

ncwriteatt(ncfile, '/', 'history', ...
    [datestr(now, 'yyyy-mm-ddTHH:MM:SSZ') ' Created by run_raw_CTD_export.m from RSK raw files']);
ncwriteatt(ncfile, '/', 'date_created',  datestr(now, 'yyyy-mm-ddTHH:MM:SSZ'));
ncwriteatt(ncfile, '/', 'date_modified', datestr(now, 'yyyy-mm-ddTHH:MM:SSZ'));

%% ════════════════════════════════════════════════════════════════════════
%  STEP 9: WRITE DATA
%  ════════════════════════════════════════════════════════════════════════
disp('  Writing data...')

ncwrite(ncfile, 'time',        time_profile);
ncwrite(ncfile, 'latitude',    ctdlat);
ncwrite(ncfile, 'longitude',   ctdlon);
ncwrite(ncfile, 'station_id',  station_char');
ncwrite(ncfile, 'n_samples',   int32(n_samp));
ncwrite(ncfile, 'sampling_frequency', freq);
ncwrite(ncfile, 'max_pressure', max_pres);
ncwrite(ncfile, 'sss',          sss);

ncwrite(ncfile, 'atm_pressure_1',   atmhPa1);
ncwrite(ncfile, 'atm_pressure_2',   atmhPa2);
ncwrite(ncfile, 'wind_speed',       atmwindspeed);
ncwrite(ncfile, 'wind_direction',   atmwinddirection);
ncwrite(ncfile, 'relative_humidity',atmrelativehumidity);

ncwrite(ncfile, 'sample_time',  t_2d);
ncwrite(ncfile, 'pressure',     P_2d);
ncwrite(ncfile, 'sea_pressure', SP_2d);
ncwrite(ncfile, 'temperature',  T_2d);
ncwrite(ncfile, 'conductivity', C_2d);

%% ════════════════════════════════════════════════════════════════════════
%  STEP 10: VERIFICATION
%  ════════════════════════════════════════════════════════════════════════
disp(' ')
disp('--- Step 10: Verification ---')

info = ncinfo(ncfile);
fprintf('  Format  : %s\n', info.Format)
fprintf('\n  Dimensions:\n')
for id = 1:length(info.Dimensions)
    fprintf('    %-14s = %d\n', info.Dimensions(id).Name, info.Dimensions(id).Length)
end
fprintf('\n  Variables (%d):\n', length(info.Variables))
for iv = 1:length(info.Variables)
    v = info.Variables(iv);
    dims_str = strjoin({v.Dimensions.Name}, ' × ');
    fprintf('    %-22s  (%s)\n', v.Name, dims_str)
end

% Roundtrip check
T_check = ncread(ncfile, 'temperature');
max_diff_T = max(abs(T_check(:) - T_2d(:)), [], 'omitnan');
fprintf('\n  Roundtrip T: max diff = %.2e\n', max_diff_T)

P_check = ncread(ncfile, 'pressure');
max_diff_P = max(abs(P_check(:) - P_2d(:)), [], 'omitnan');
fprintf('  Roundtrip P: max diff = %.2e\n', max_diff_P)

sid = ncread(ncfile, 'station_id');
fprintf('  First IDs: ');
for k = 1:min(5, np), fprintf('%s  ', strtrim(sid(:,k)')); end
fprintf('\n')
fprintf('  Last IDs: ');
for k = np-2:np, fprintf('%s  ', strtrim(sid(:,k)')); end
fprintf('\n')

n_samp_check = ncread(ncfile, 'n_samples');
fprintf('  n_samples min/max: %d / %d\n', min(n_samp_check), max(n_samp_check))

sss_check = ncread(ncfile, 'sss');
n_sss_valid = sum(~isnan(sss_check));
fprintf('  SSS valid: %d / %d profiles\n', n_sss_valid, np)

d = dir(ncfile);
fprintf('\n  File size: %.1f KB (%.1f MB)\n', d.bytes/1024, d.bytes/1024/1024)

disp(' ')
disp('╔════════════════════════════════════════════════════════════════════╗')
disp('║  EXPORT RAW NETCDF — DONE                                          ║')
disp('╚════════════════════════════════════════════════════════════════════╝')
fprintf('  File: %s\n', ncfile)
fprintf('  %d profils | obs_max = %d | 5 variables (obs×profile)\n', np, obs_max)
disp(' ')
disp('  For inspection:')
fprintf('    ncdisp(''%s'')\n', ncfile)
disp('  For CF validation:')
fprintf('    cfchecks %s\n', ncfile)
