% GNSS-DPE tool main script %
% Adri Gusi, Pau Closas, Shuo Tang
% v3.0 last update 04/18/2022
% This script generates GPS signals and computes the receiver position with
% DPE based algorithms with Doppler shift and carrier phase.
% Current implementation:
% -7 GPS satellites at a fixed location
% -x(t) = As(t)*D*C 
% -Non-coherent and coherent integration methods
% -DPE algorithm with ARS

close all;
clear all; clc; format long;

% tic
%% load configuration file
config = 'ConfigFile';
eval(config)
% load true parameters
trueParam.UserPosition = param.UserPosition;
trueParam.UserVelocity = param.UserVelocity;
trueParam.deltaT = param.deltaT;
trueParam.deltaTdot = param.deltaTdot;
%% number of samples calculation
% Number of samples of the Local Replica
NsamplesLocal=CodePeriod*fs*CoherentIntegrations;   
 %Number of samples of the Received Signal (Data)
NsamplesData=CodePeriod*fs*CoherentIntegrations*NonCoherentIntegrations;   

%% Memory allocation.
PosErrLS = zeros(length(CNosim),Nexpe);
PosErrDPE = zeros(length(CNosim),Nexpe);
CBErrDPE = zeros(length(CNosim),Nexpe);
PosEstDPE = zeros(Nexpe,3,length(CNosim));

%% Signal generation (sigen struct)
fprintf('Generating signals...')
caCode_filt = filtcaCode(config);
sigen = signalGen0(config,caCode_filt,trueParam);

fprintf('Done\n')

%% Start Simulation
for CNo_idx=1:length(CNosim)
%     tic
   parfor exp_idx=1:Nexpe
        fprintf('Simulation for CNo = %d, Experiment #%d \n',CNosim(CNo_idx),exp_idx)
        CNo=CNosim(CNo_idx)*ones(numSV,1);
        %% Signal + noise
        x_delay_noise = receivedSignal0(sigen,config,CNo);

        %% DPE approach ARS
        [PosErrDPE(CNo_idx,exp_idx),CBErrDPE(CNo_idx,exp_idx),PosEstDPE(exp_idx,:,CNo_idx)] = ...
            DPEarsPVT(x_delay_noise,config,caCode_filt,CNo_idx);
        fprintf('Done\n');
    end
%     toc
end
% toc
% Compute RMSEs
RMSE_DPE = sqrt(mean(PosErrDPE.^2,2));
RMSECB_DPE = sqrt(mean(CBErrDPE.^2,2));
%% Ziv-Zakai bound computation
[fZZLB_DPE, fZZLB_CB] = computeZZBDPE_4d(config,sigen);
% fZZLB_DPE = computeZZBDPE(config,sigen);

%% PLOTS
figure,h=semilogy(CNosim0,fZZLB_DPE,'r-');
legend('ZZB DPE')
grid
set(h,'Linewidth',2)
figure,h=semilogy(CNosim0,fZZLB_CB,'r-');
legend('ZZB DPE')
grid
set(h,'Linewidth',2)
% figure,h=semilogy(CNosim,RMSE_DPE,'o-');
% legend('MLE DPE')
% grid
% set(h,'Linewidth',2)
if DPEflag == 0
    figure,
    h=semilogy(CNosim,RMSE_LS,'b-',CNosim0,fZZLB_2SP,'r-');
    legend('MLE 2SP','ZZB 2SP')
    grid
    set(h,'Linewidth',2)
else 
    figure,
    h=semilogy(CNosim,RMSE_DPE,'b-',CNosim0,fZZLB_DPE,'r-');
    legend('MLE DPE','ZZB DPE')
    grid
    set(h,'Linewidth',2)
    xlabel('CN0 [dB-Hz]')
    ylabel('RMSE [m]')
    figure,
    h1=semilogy(CNosim,RMSECB_DPE,'b-',CNosim0,fZZLB_CB,'r-');
    legend('MLE DPE','ZZB DPE')
    grid
    set(h1,'Linewidth',2)
    xlabel('CN0 [dB-Hz]')
    ylabel('RMSE [s]')
end
