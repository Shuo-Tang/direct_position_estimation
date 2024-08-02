%% Configuration file for parameters adjustment


%% Constants.
c = 299792458;                                  %   Speed of light.
f0 = 1575.42e6;                                  %   Carrier frequency (L1-E1 band).
lambda = c/f0;
%% $Signal Parameters
% GPS E1
CodePeriod  =   1e-3;               % Code Period
CoherentIntegrations= 1;            % Number of Coherent Integrations
NonCoherentIntegrations = 1;        % Number of Non Coherent Integrations
m=1;                                % m
n=1;                                % n
type='BPSK';                        % BPSK, BOCcos, BOCsin
Tc=1/(n*1.023e6);                   % Chip time
Ts=1/(2*m*1.023e6);
fs=50e6;                            % Sampling frequency
dt=1/fs;                            % Sampling period
CodeLen = 1023;
Tchip = CodePeriod / CodeLen;
%% Scenario Parameters
% Set GPS SVs&user's position and other parqameters
load('param.mat')
numSV = 7;
SatPRN=  [12    15    17    19    24    25    32];
DPEflag = 1;
DPEcase = 2;% 0- no Doppler, no Carrier
            % 1- Doppler, no Carrier
            % 2- Doppler and Carrier

%% Simulation parameters
% Set C/N0 (dB Hz) of each SV signal
CNosim0 = (30:1:60);
CNosim = (30:1:60);
% Number of experiments for each simulated C/N0
Nexpe = 100;
% Simulate MLE flag
simulate_mle = 1;
% Plot bound flag
compute_zzb = 1;


%% 2-steps parameters
% LS iterations
num2stepsIterations = 10;

%% DPE parameters (ARS)

Niter=1800;
gamma_est=zeros(3,Niter+1);
EstRxClkBias=zeros(1,Niter+1);
% EstRxClkBias=dt*1.3;
contraction = 2;          % contraction parameter

% Dmax0 = [1e2*ones(1,2),1e1*ones(1,3),1e-2*ones(1,2),1e-3*ones(1,8),1e-3*ones(1,10),5e-4*ones(1,6)]; % v1.0  
% Dmin0 = [1e1*ones(1,2),1*ones(1,3),1e-3*ones(1,2),1e-4*ones(1,8),1e-4*ones(1,10),5e-5*ones(1,6)]; % v1.0 
% Dmax0 = [5e2*ones(1,2),1e2*ones(1,3),1e-2*ones(1,4),5e-3*ones(1,2),1e-3*ones(1,4),1e-3*ones(1,10),5e-4*ones(1,6)]; % v1.1 
% Dmin0 = [5e1*ones(1,2),5*ones(1,3),1e-3*ones(1,4),5e-4*ones(1,2),1e-4*ones(1,4),1e-4*ones(1,10),5e-5*ones(1,6)]; % v1.1
Dmax0 = [1e2*ones(1,2),1e0*ones(1,2),1e-1,5e-2,1e-2*ones(1,5),1e-3*ones(1,4),1e-3*ones(1,10),5e-4*ones(1,6)]; % v1.2 
Dmin0 = [1e1*ones(1,2),1e-1*ones(1,2),1e-2,5e-3,1e-3*ones(1,5),1e-4*ones(1,4),1e-4*ones(1,10),5e-5*ones(1,6)]; % v1.2

Dmax = Dmax0(find(CNosim0==CNosim(1)):find(CNosim0==CNosim(end)));
Dmin = Dmin0(find(CNosim0==CNosim(1)):find(CNosim0==CNosim(end)));

% Dmax_clk0 = [1e1*ones(1,2),1e0*ones(1,2),1e-2,1e-3,1e-5,1e-6,1e-7,1e-9,1e-10,1e-12*ones(1,20)]; % v1.0 
% Dmin_clk0 = [1e0*ones(1,2),1e-1,1e-2,1e-3,1e-4,1e-6,1e-7,1e-8,1e-10,1e-11,1e-13*ones(1,10),1e-14*ones(1,10)]; % v1.0 
Dmax_clk0 = [1e1*ones(1,2),1e1*ones(1,2),5e-1,5e-7,5e-9,1e-10*ones(1,4),1e-12*ones(1,20)]; % v1.1 
Dmin_clk0 = [1e0*ones(1,2),1e-1*ones(1,2),5e-2,5e-8,5e-10,1e-11*ones(1,4),1e-13*ones(1,10),1e-14*ones(1,10)]; % v1.1
Dmax_clk = Dmax_clk0(find(CNosim0==CNosim(1)):find(CNosim0==CNosim(end)));
Dmin_clk = Dmin_clk0(find(CNosim0==CNosim(1)):find(CNosim0==CNosim(end)));

if length(CNosim) ~= length(Dmax)
    error('The Dmax setting size does not match the size of CNo.')
end
% dmax = 0.01;
% dmin = 0.0001;
% dmax_clk=dt/10;
% dmin_clk=dt/1000;
%% DPE parameters (GS)
RangeSet = 25*ones(1,31);
switch DPEcase
    case 0
        DensitySet = [40*ones(1,2),15,5,2*ones(1,12),1*ones(1,15)];
    case 1
        DensitySet = [40*ones(1,2),10,5,2*ones(1,12),1*ones(1,15)];
    case 2
        DensitySet = [40*ones(1,2),10,4,1,0.04,0.004,0.001*ones(1,8),4e-4*ones(1,16)];
end




