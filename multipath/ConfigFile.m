%% Configuration file for parameters adjustment


%% Constants.
c = 299792458;                                  %   Speed of light.
f0 = 1575.42e6;                                  %   Carrier frequency (L1-E1 band).
lambda = c/f0;
%% $Signal Parameters
% GPS L1
CodePeriod  =   1e-3;               % Code Period
CoherentIntegrations= 1;            % Number of Coherent Integrations
NonCoherentIntegrations = 1;        % Number of Non Coherent Integrations
m=1;                                % m
n=1;                                % n
type='BPSK';                        % BPSK, BOCcos, BOCsin
Tc=1/(n*1.023e6);                   % Chip time
Ts=1/(2*m*1.023e6);
fs=25375e3;                            % Sampling frequency
dt=1/fs;                            % Sampling period
travelTimeOffset = 68.802*1e-3;
Tchip = CodePeriod /1023;
% IF = 420e3;
%% Scenario Parameters
% Set GPS SVs&user's position and other parqameters
% load('param.mat')
numSV = 7;
DPEflag = 1;

%% Simulation parameters
% Set C/N0 (dB Hz) of each SV signal
CNosim = 50;
% Number of experiments for each simulated C/N0
Nexpe = 150;
% Simulate MLE flag
simulate_mle = 1;
% Plot bound flag
compute_zzb = 1;


%% DPE parameters (ARS)

Niter=5000;
gamma_est=zeros(3,Niter+1);
EstRxClkBias=zeros(1,Niter+1);
% EstRxClkBias=dt*1.3;
contraction = 2;          % contraction parameter
Dmax = 20;
Dmin = 0.1;
Dmax_clk=dt/100;
Dmin_clk=dt/1000;
%% DPE parameters (GS)
% RangeSet = 25*ones(1,31);
% DensitySet = [40*ones(1,2),10,4,1,0.04,0.004,0.001*ones(1,8),4e-4*ones(1,16)];

SearchRange = 100;
SearchDensity = 2;

posSearchNum = 100;
posSearchDensity = 1;
cbSearchNum = 20;
cbSearchDensity = 1e-8;



