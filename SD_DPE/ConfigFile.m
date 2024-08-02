%% Configuration file for parameters adjustment


%% Constants.
c = 299792458;                                  %   Speed of light.
f0 = 1575.42e6;                                  %   Carrier frequency (L1-E1 band).
lambda = c/f0;
%% Signal Parameters
% GPS E1
CodePeriod  =   1e-3;               % Code Period
CoherentIntegrations= 1;            % Number of Coherent Integrations
NonCoherentIntegrations = 1;        % Number of Non Coherent Integrations
m=1;                                % m
n=1;                                % n
type='BPSK';                        % BPSK, BOCcos, BOCsin
Tc=1/(n*1.023e6);                   % Chip time
Ts=1/(2*m*1.023e6);
fs=40.92e6;                            % Sampling frequency
dt=1/fs;                            % Sampling period
CodeLen = 1023;
Tchip = CodePeriod / CodeLen;
%% Scenario Parameters
% Set GPS SVs&user's position and other parqameters
% load('param.mat')
% numSV = 7;
% SatPRN=  [12    15    17    19    24    25    32];
% SatPRN=  [12    15    18    19    24    28    32];
% DPEflag = 1;
% DPEcase = 2;% 0- no Doppler, no Carrier
%             % 1- Doppler, no Carrier
%             % 2- Doppler and Carrier

%% Simulation parameters




% Set C/N0 (dB Hz) of each SV signal
% CNosim0 = (30:1:60);
% CNosim = 50;
% % Number of experiments for each simulated C/N0
% Nexpe = 100;
% % Simulate MLE flag
% simulate_mle = 1;
% % Plot bound flag
% compute_zzb = 1;


%% 2-steps parameters
% LS iterations
num2stepsIterations = 10;

%% DPE parameters (ARS)





