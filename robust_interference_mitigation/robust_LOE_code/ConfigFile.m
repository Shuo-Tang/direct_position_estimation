%% Configuration file for parameters adjustment


%% Constants.
c                                                   =   299792458;                                  %   Speed of light.
f0                                                  =   1575.42e6;                                  %   Carrier frequency (L1-E1 band).

%% $Signal Parameters
% GPS E1
CodePeriod  =   1e-3;               % Code Period
CoherentIntegrations= 1;            % Number of Coherent Integrations
NonCoherentIntegrations = 1;        % Number o Non Coherent Integrations
m=1;                                % m
n=1;                                % n
type='BPSK';                        % BPSK, BOCcos, BOCsin
Tc=1/(n*1.023e6);                   % Chip time
Ts=1/(2*m*1.023e6);
fs=50e6;                            % Sampling frequency
dt=1/fs;                            % Sampling period
fn=2e6;
order=36;
%% Scenario Parameters
% Set GPS SVs and user positions
% numSV=7;
% load('SatPositions.mat')
% SatPosition=corrSatPosition(1:7, :);
% UserPosition=[3.915394273911475e+06 2.939638207807819e+05 5.009529661006817e+06];
% SatPRN=  [12    15    17    19    24    25    32 ];
load('param.mat')
numSV=size(param.SatPosition,1); 
UserPosition = param.UserPosition;
SatPosition = param.SatPosition;
SatPRN = 1:numSV;
%% Simulation parameters
% Set C/N0 (dB Hz) of each SV signal
CNosim = 44;%45:1:50;
% Number of experiments for each simulated C/N0
Nexpe = 5e4;%5e4;%2e3;
% Simulate MLE flag
simulate_mle = 1;
% Plot bound flag
compute_zzb = 1;
% Apply RIM flag
RIMuse = 1;
% EstimateTrueNoise
estimateTrueNoise = 1;
% Plot estimated CN0
plot_estimated_cn0 = 1;

%% 2-steps parameters
% LS iterations
num2stepsIterations = 10;

%% DPE parameters (ARS)

Niter=3e4;%2e3;%1e5;
gamma_est=zeros(3,Niter+1);
amp_est = zeros(numSV, Niter+1);
EstRxClkBias=zeros(1,Niter+1);
% EstRxClkBias=dt*1.3;
contraction = 2;          % contraction parameter
dmax = 10;%1e4;
dmin = 1e-1;%0.01;
dmax_clk=dt/10;
dmin_clk=dt/100;


NormalizaFactor = sqrt(NonCoherentIntegrations)*CoherentIntegrations*CodePeriod*fs;
CN0_est_ind = 1;

