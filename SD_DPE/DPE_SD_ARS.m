clear
clc
close all

%% Set parameters and positions
config = 'ConfigFile';
eval(config)
nSamples_1ms = fs * 1e-3;
plot_delays = 0;
noise_on = 0;
iono_tropo_on = 0;
range_error_on = 0;
plot_CAFs = 0;
DPE_SD_enable = 1;
DPE_standard_enable = 1;
% === Reference station =====
ref_info = zeros(1, 8);
ref_info(1:3) = [4777973.177,176346.307,4207663.620]; % position
ref_info(4) = 0; % clock bias
ref_info(5:7) = [0, 0, 0]; % velocity
ref_info(8) = 0; % clock drift

% === Rover =====
baseline = [1000, 1000, 800];
rov_ecef = ref_info(1:3) + baseline;
rov_info = zeros(1, 8);
rov_info(1:3) = rov_ecef; % position
rov_info(4) = 0; % clock bias
rov_info(5:7) = [0, 0, 0]; % velocity
rov_info(8) = 0; % clock drift

% === Satellite =====
sat_info_original = load("load_data/sat8_info_multipath_fix_pos.mat");
sat_info_original = sat_info_original.sat_info(:,:,1);
nSat = size(sat_info_original, 1);
sat_cb = zeros(nSat, 1);
sat_cd = zeros(nSat, 1);
sat_info = [sat_info_original(:,1:3), sat_cb,...
    sat_info_original(:,4:6), sat_cd, sat_info_original(:,9)];
caCode = zeros(32, 1023);
for i = 1:32
    caCode(i,:) = genCAcode(i);
end
%% Generate two signals and observations
iono_error = 10;
tropo_error = 3;
iono_sigma = 1;
tropo_sigma = 0.1;
CNo = 35;
% === Reference station =====
if iono_tropo_on == 1
    iono_ref = iono_error + iono_sigma*randn(nSat, 1);
    tropo_ref = tropo_error + tropo_sigma*randn(nSat, 1);
    % iono_ref = load("load_data/iono_ref.mat");
    % iono_ref = iono_ref.iono_ref;
    % tropo_ref = load("load_data/tropo_ref.mat");
    % tropo_ref = tropo_ref.tropo_ref;
else
    iono_ref = zeros(nSat, 1);
    tropo_ref = zeros(nSat, 1);
end
[x_ref, pseudorange_true_ref] =...
    signalGenerate_iono(config, ref_info, sat_info, caCode, iono_ref, tropo_ref);

% === Rover station ===== 
if iono_tropo_on == 1
    iono_rov = iono_error + iono_sigma*randn(nSat, 1);
    tropo_rov = tropo_error + tropo_sigma*randn(nSat, 1);
    % iono_rov = load("load_data/iono_rov.mat");
    % iono_rov = iono_rov.iono_rov;
    % tropo_rov = load("load_data/tropo_rov.mat");
    % tropo_rov = tropo_rov.tropo_rov;
else
    iono_rov = zeros(nSat, 1);
    tropo_rov = zeros(nSat, 1);
end
[x_rov, pseudorange_true_rov] = signalGenerate_iono(config, rov_info, sat_info, caCode, iono_rov, tropo_rov);
if noise_on == 1
    raw_signal_rov = receivedSignal(x_rov, config, CNo);
else
    raw_signal_rov = x_rov;
end

%% Correlate between raw signal and reconstructed signal
if range_error_on == 1
    pseudorange_est_ref = pseudorange_true_ref + 1*randn(nSat, 1);
else
    pseudorange_est_ref = pseudorange_true_ref;
end


%% implement DPE-SD -- grid search
if DPE_SD_enable == 1
% === settings =====
% set ARS parameter
dMax = 20;
dMin = 0.01;
contraction = 2;
nIter = 400;
b_initial = [baseline(1:2) + 50*randn(1,2), baseline(3)];
% === implement DPE-SD =====
[b_est, ~] = arsSearch_baseline(config, raw_signal_rov, b_initial, ref_info, sat_info, caCode,...
        pseudorange_est_ref, nIter, dMax, dMin, contraction); 
% === record sol =====
sol_sd = b_est

end

%% implement standard DPE -- grid search
if DPE_standard_enable == 1
% set ARS parameter
dMax = 20;
dMin = 0.01;
contraction = 2;
nIter = 400;
pos_initial = [rov_info(1:2) + 50*randn(1,2), rov_info(3)];
rov_initial = [pos_initial(1:3), rov_info(4:end)];
% === implement standard DPE =====
[pos_est, ~] = arsSearch(config, raw_signal_rov, rov_initial, sat_info, caCode,...
        nIter, dMax, dMin, contraction); 
% === look for maximum =====
sol_pos = pos_est
end