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

rng(100)
%% load configuration file
config = 'ConfigFile';
eval(config)

G_rx=10^(2.5); %G_rx= 25 dB
G_tx= 10^(0); %G_tx= 0 dB
c=3*10^8;
lambda=c/1170e6;
Distance=[logspace(1,1.5,10),logspace(1.6,4,10)];
L_Friis=G_rx*G_tx.*(lambda/(4*pi)./Distance).^2;
P_tx=10^(106/10)/126; %P_tx= 106 dBW with P_N assumed as 126, keep same as DD simulation.
P_N=1; 
JN_Array=10*log10(P_tx/P_N.*L_Friis);
JN_Array=[-20:5:0, 3:3:35, 36:3:50];
% JN_Array=[-20:5:0, 3:3:30];
%% Memory allocation.
PosErrLS=zeros(length(JN_Array),Nexpe);
PosErrDPE=zeros(length(JN_Array),Nexpe);
Pos_est=zeros(length(JN_Array),3, Nexpe);
cn0= zeros(length(JN_Array),numSV,Nexpe);
%% Signal generation (sigen struct)
sigen = signalGen(config);
% meanNoise = computeMeanNoise(config,sigen);

if simulate_mle
    %% Start Simulation
    for JN_idx=1:length(JN_Array)
%         Threshold(CNo_idx)
         JN = JN_Array(JN_idx);
         fprintf('Simulation for JN_idx = %d \n',JN_idx)
        parfor exp_idx=1:Nexpe
            CNo=CNosim*ones(numSV,1);
            
            %% Signal + noise 
            x = receivedSignal(sigen,config,CNo);
            x_CW = Jam_Data_Generator(config);
            x = x + x_CW.*sqrt(10^(JN/10)*P_N);
            
%             x_DME = DME_generator(config);
%             x = x + x_DME.*sqrt(10^(JN/10)*P_N);
            %% Apply RIM
            if RIMuse
%                 x = RIM_TD(x, 1.345);
                x = RIM_FD(x, 1.345);
            end
            
            %% Perform coherent/non-coherent integration times
            r = correlateSignal(sigen,config,x);
            
            %% Estimate CN0
%             cn0(CNo_idx,:,exp_idx) = estimateCn0(r,config,meanNoise);
            
            %% 2-steps: Conventional approach estimation
%             PosErrLS(JN_idx,exp_idx) = conv2stepsPVT(r,config);
                                               
            %% DPE approach ARS (accelerated random search)
            [PosErrDPE(JN_idx,exp_idx), Pos_est(JN_idx,:,exp_idx)] = DPEarsPVT(r,config);
%             if rem(exp_idx, 1000) == 0
%                 fprintf('Simulation for JN_idx = %d, exp_idx = %d \n',JN_idx, exp_idx)
%             end
        end
    end
    
    % Compute RMSEs
    RMSE_LS=sqrt(mean(PosErrLS.^2,2));
    RMSE_DPE=sqrt(mean(PosErrDPE.^2,2));
    averageCn0= (mean(mean(cn0,2),3));
    save('PosErrDPE_CW_FDRIM.mat', 'CNosim', 'JN_Array', 'PosErrDPE', 'RMSE_DPE');
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
