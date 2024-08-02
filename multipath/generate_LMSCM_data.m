clear
close all
clc;
clear classes

%% Load Predefined Setting
pth = cd;
addpath ([pth,'\Knife Edge']);
addpath ([pth,'\Parameter Files']);

%% Load Settings
% predefined parameter
% user parameters
X = 4777973.177; Y = 176346.307; Z = 4207663.62;
user_pos_ecef_initial = [X, Y, Z];
user_pos_lla_initial = ecef2lla(user_pos_ecef_initial);
user_clock_bias = 2e-8;
user_pos = user_pos_ecef_initial;
c = 299792458;
% receiver parameters
fs = 20.46e6;      % sampling frequency
ts = 1/fs;
duration = 10;  % sec
samplesPerCode = fs * 1e-3;
tc = 1/1.023e6; % period of one chip 
codeLength = 1023;
fb = 50;
tb = 1 / fb;
CN0 = 40;
snr = 10^((CN0-10*log10(fs))/10);
noise_power = 1;
signal_power = noise_power * snr;
alpha = sqrt(signal_power / 1);
% satellite parameters
sat_info = load("load_data/sat8_info_multipath.mat");
sat_info = sat_info.sat_info;
prn = sat_info(:,9,1);
nSat = length(prn);

%% Generate LMSCM Delays
for iSat = 1:nSat
    iSat
    validate_LOS = false;
    % === LMSCM settings =====
    Parameters.SampFreq = 100;        % Hz
    Parameters.MaximumSpeed = 5;     % km/h
    Parameters.NumberOfSteps = duration * Parameters.SampFreq;

    user_type = 'Pedestrian';
    scenario = 'Suburban';

    Parameters.SatAzimut = sat_info(iSat, 7, 1);       % Deg
    Parameters.SatElevation = sat_info(iSat, 8, 1);
    % === Generate delays vecotrs for each satellite =====
    while validate_LOS == false
        [LOS_amplitudes, LOS_delays, echo_amplitudes, echo_delays] =...
            generate_delays_and_CIRs(user_type, scenario, Parameters);
        try
            LOS_amp_array = abs(cell2mat(LOS_amplitudes)); 
            validate_LOS = all(LOS_amp_array(:) > 0.8);
        catch
            disp("more than one LOS\n")
        end
    end
    % [LOS_amplitudes, LOS_delays, echo_amplitudes, echo_delays] =...
    %         generate_delays_and_CIRs(user_type, scenario, Parameters);
    % for i = 1: Parameters.NumberOfSteps
    %     LOS_amplitudes{i} = 0.01*randn;
    % end
    % === test maximum delay =====
    max_delay = 0;
    for i = 1:Parameters.NumberOfSteps
        max_delay = max([max_delay, max(echo_delays{i})]);
    end
    lengthDelayVec = round(max_delay * fs) + 1;
    % === convert to delay sampleing bins =====
    delay_loss = zeros(Parameters.NumberOfSteps, lengthDelayVec);
    for i = 1:Parameters.NumberOfSteps
        % LOS amplitude
        delay_loss(i, 1) = sum(LOS_amplitudes{i});
        % NLOS amplitude
        delays = echo_delays{i};
        delay_amplitudes = echo_amplitudes{i};
        nDelay = length(delays);
        for iDelay = 1: nDelay
            delay_ind = round(delays(iDelay) / ts);
            delay_loss(i, delay_ind + 1) = delay_loss(i, delay_ind + 1) + delay_amplitudes(iDelay);
        end
    end
    save_path = sprintf("load_data/delay_loss_LOS/delay_loss_sat%d_%s_%s.mat",...
        prn(iSat), user_type, scenario);
    save(save_path, "delay_loss")
end



%% NLOS
% user_type = 'Car';
% scenario = 'Urban';
% for iSat = 1:3
%     load_path_NLOS = sprintf("load_data/delay_loss_3NLOS/delay_loss_sat%d_%s_%s_NLOS.mat",...
%     prn(iSat), user_type, scenario);
%     load_path_LOS = sprintf("load_data/delay_loss_3NLOS/delay_loss_sat%d_%s_%s.mat",...
%     prn(iSat), user_type, scenario);
%     delay_loss_NLOS = load(load_path_NLOS);
%     delay_loss_LOS = load(load_path_LOS);
%     delay_loss_NLOS = delay_loss_NLOS.delay_loss;
%     delay_loss_LOS = delay_loss_LOS.delay_loss;
% 
%     delay_loss = delay_loss_LOS;
%     n_LOS = size(delay_loss_LOS, 2);
%     n_NLOS = size(delay_loss_NLOS, 2);
%     n = min(n_LOS, n_NLOS);
% 
%     delay_loss(51:450, 1:n) = delay_loss_NLOS(51:450, 1:n);
%     save_path = sprintf("load_data/delay_loss_3NLOS/delay_loss_sat%d_%s_%s.mat",...
%          prn(iSat), user_type, scenario);
%     save(save_path, "delay_loss")
% 
% 
% end