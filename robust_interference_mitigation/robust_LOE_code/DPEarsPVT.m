function [PosErrDPE, Pos_est] = DPEarsPVT(r,config)

%% load configuration file
eval(config)

randomDelay = 0; %%% magic number....
gamma_est = zeros(3, Niter);
%% memory allocation
EstRange=zeros(1,numSV);
MaxCorr=zeros(1,numSV);

gamma_est(:,1) = UserPosition+100*(2*rand(3,1)-1)';
% gamma_est(3, 1) = UserPosition(1,3);
EstRxClkBias(:,1)=-randomDelay*Tc;


for kSV                                         =   1 : numSV
    % Compute range estimate with corrected satellite position, atmosphere
    % corrections and satellite clock error.
    EstRange(kSV)                               =   norm(SatPosition(kSV,:) - gamma_est(:,1)');
end;

EstFracDelay=mod(EstRange/c+EstRxClkBias(:,1)+dt,1e-3);
% EstFracDelay=EstFracDelay+EstRxClkBias(:,1);
% EstFracDelay(EstFracDelay<0)=EstFracDelay(EstFracDelay<0)+1e-3;
% EstFracDelay(EstFracDelay>1e-3)=EstFracDelay(EstFracDelay>1e-3)-1e-3;

for kSV                                         =   1 : numSV
    
    %     round(FracDelay*Fs)
    aux=round(EstFracDelay(kSV)*fs);
    aux(aux==0)=1e-3*fs;
    MaxCorr(kSV)=r(kSV,aux);
    
end
J_ant=sum(MaxCorr);
amp_est(:,1) = MaxCorr./NormalizaFactor^2;

d = dmax;
d_clk=dmax_clk;

for it = 1:Niter-1        %%% ARS algorithm iterations
    
    % draw a random movement
    rand_point = gamma_est(:,it) + d*(2*rand(3,1)-1);
%     rand_point(3,1) = UserPosition(1,3);
    % rand_clk = EstRxClkBias(:,it)+d_clk*(2*rand-1);
    rand_clk = 0;
    for kSV                                         =   1 : numSV
        
        EstRange(kSV)                               =   norm(SatPosition(kSV,:) - rand_point');
    end
    
    EstFracDelay=mod(EstRange/c+rand_clk+dt,1e-3);
    
    for kSV                                         =   1 : numSV
        
        aux=round(EstFracDelay(kSV)*fs);
        aux(aux==0)=1e-3*fs;
        MaxCorr(kSV)=r(kSV,aux);
        
    end
    
    
    J = sum(MaxCorr);
    
    % select or discard point
    if J > J_ant
        gamma_est(:,it+1) = rand_point;
        EstRxClkBias(:,it+1)=rand_clk;
        amp_est(:, it+1) = MaxCorr./NormalizaFactor^2;
        J_ant = J;
        d = dmax;
        d_clk=dmax_clk;
    else
        gamma_est(:,it+1) = gamma_est(:,it);
        EstRxClkBias(:,it+1)=EstRxClkBias(:,it);
        amp_est(:,it+1) = amp_est(:,it);
        d = d/contraction;
        d_clk=d_clk/contraction;
    end
    if d < dmin
        d = dmax;
    end
    
    if d_clk < dmin_clk
        d_clk =dmax_clk;
    end
end
% DPE position estimation
PosErrDPE=norm(gamma_est(:,it+1)'-UserPosition);
Pos_est = gamma_est(:,it+1);
% CN0_est = 10.*log10(amp_est(:, it+1).*(2*fn)); 
% CN0_est = CN0_est (CN0_est_ind);