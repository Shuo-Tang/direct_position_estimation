% GNSS-DPE tool main script %
% Adrià Gusi, Pau Closas
% v2.0
% This script generates GPS signals and computes the receiver position with
% two-steps and DPE based algorithms.
% Current implementation:
% -7 GPS satellites at a fixed location
% -Code based only
% -Non-coherent and coherent integration methods
% -LS method implemented with Integer-Millisecon Rollover correction
% -DPE algorithm with ARS implementation

% close all;
clear all; 
% clc; 
format long;


%% load configuration file
config = 'ConfigFile';
eval(config)

%% Memory allocation.
PosErrLS=zeros(length(CNosim),Nexpe);
PosErrDPE=zeros(length(CNosim),Nexpe);
CN0_est=zeros(length(CNosim),Nexpe);
cn0= zeros(length(CNosim),numSV,Nexpe);
%% Signal generation (sigen struct)
sigen = signalGen(config);
meanNoise = computeMeanNoise(config,sigen);

if simulate_mle
    %% Start Simulation
    for CNo_idx=1:length(CNosim)
        CNosim(CNo_idx)
        
        for exp_idx=1:Nexpe
            CNo=CNosim(CNo_idx)*ones(numSV,1);
            
            %% Signal + noise 
            x = receivedSignal(sigen,config,CNo);
            
            
            %% Apply RIM
            if RIMuse
                x = RIM(x);
            end
            
            %% Perform coherent/non-coherent integration times
            r = correlateSignal(sigen,config,x);
            
            %% Estimate CN0
            cn0(CNo_idx,:,exp_idx) = estimateCn0(r,config,meanNoise);
            
            %% 2-steps: Conventional approach estimation
            PosErrLS(CNo_idx,exp_idx) = conv2stepsPVT(r,config);
                                               
            %% DPE approach ARS (accelerated random search)
            [PosErrDPE(CNo_idx,exp_idx), CN0_est(CNo_idx,exp_idx)] = DPEarsPVT(r,config);

        end
    end
    
    % Compute RMSEs
    RMSE_LS=sqrt(mean(PosErrLS.^2,2));
    RMSE_DPE=sqrt(mean(PosErrDPE.^2,2));
    averageCn0= (mean(mean(cn0,2),3));
    
end

if compute_zzb
    %% Ziv-Zakai bound computation
    computeZZB
end

%% PLOTS
if simulate_mle || compute_zzb
    figure,
    if simulate_mle && compute_zzb
        h=semilogy(CNosim,RMSE_LS,'b-.',CNosim,RMSE_DPE,'b',CNosim,fZZLB_2SP,'r-.',CNosim,fZZLB_DPE,'r');
        legend('MLE 2SP','MLE DPE','ZZB 2SP','ZZB DPE')
        grid
        set(h,'Linewidth',2)
    elseif simulate_mle
        h=semilogy(CNosim,RMSE_LS,'b-.',CNosim,RMSE_DPE,'b');
        legend('MLE 2SP','MLE DPE')
        grid
        set(h,'Linewidth',2)
    else
        h=semilogy(CNosim,fZZLB_2SP,'r-.',CNosim,fZZLB_DPE,'r');
        legend('ZZB 2SP','ZZB DPE')
        grid
        set(h,'Linewidth',2)
    end
    xlabel('CN0 [dB-Hz]')
    ylabel('RMSE [m]')
end


if plot_estimated_cn0
    figure
    h=plot(CNosim,averageCn0,CNosim,CNosim);
    legend('Estimated CN0','True CN0')
    grid
    set(h,'Linewidth',2)
end
