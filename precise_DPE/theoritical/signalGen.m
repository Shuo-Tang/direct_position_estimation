function [sigen] = signalGen(config,caCode_filt,gamma)


%% load configuration file
eval(config)
dimen = 3;
% fixed parameters
SatPosition = param.SatPosition(:,1:dimen);
SatVelocity = param.SatVelocity(:,1:dimen);
deltat = param.deltat;
% estimated parameters
UserPosition = gamma.UserPosition(:,1:dimen);
UserVelocity = gamma.UserVelocity(:,1:dimen);
deltaT = gamma.deltaT;
deltaTdot = gamma.deltaTdot;
UnitDirection = (UserPosition-SatPosition)./vecnorm(UserPosition-SatPosition,2,2);

%% number of samples calculation
NsamplesLocal=CodePeriod*fs*CoherentIntegrations;   % Number of samples of the Local Replica
NsamplesData=CodePeriod*fs*CoherentIntegrations*NonCoherentIntegrations;    %Number of samples of the Received Signal (Data).


%% memory allocation
GeoRange = zeros(1,numSV);
Range = zeros(1,numSV);
FracDelay = zeros(1,numSV);
fd = zeros(1,numSV);
phi = zeros(1,numSV);
%% Compute range and fractional delays for each SV
for kSV=1:numSV 
    GeoRange(kSV) = norm(SatPosition(kSV,:) - UserPosition); 
    Range(kSV) = GeoRange(kSV) + c*(deltaT - deltat(kSV));  
    FracDelay(kSV) = mod(Range(kSV)/c,CodePeriod);
    fd(kSV) = -(UserVelocity-SatVelocity(kSV,:))*UnitDirection(kSV,:)'*f0/c...
        *(1+deltaTdot);
    phi(kSV) = 2*pi*(Range(kSV)/lambda);
end
%% Signal Generated according to gamma
    %% memory set
    x_delay1 = zeros(numSV,NsamplesData);
    x_delay2 = zeros(numSV,NsamplesData);
    x_delay = zeros(numSV,NsamplesData);
    s = zeros(numSV,NsamplesData);
    PrevNCOIndex = -FracDelay/Tc;
    % setting for filter
%     fn=2e6;
%     order=9;
%     wn=pi*fn/(pi*fs/2);
%     h=fir1(order,wn);
    %% generate signal
    for kSV = 1:numSV
%         Code = genCAcode(SatPRN(kSV));
%         ii = 1 : NsamplesData;
%         %% generate C/A code, Doppler and carrier phase
%         switch DPEcase
%             case 0
%                 x_delay1(kSV,:) = Code((1 + mod(round(PrevNCOIndex(kSV)+(ii - 1)/fs/Tchip), CodeLen)));
%             case 1
%                 x_delay1(kSV,:) = Code((1 + mod(round(PrevNCOIndex(kSV)+(ii - 1)/fs/Tchip), CodeLen))).*...
%                     exp(2*pi*fd(kSV)*(1:1:NsamplesData)*dt*1i);
%             case 2
%                 x_delay1(kSV,:) = Code((1 + mod(round(PrevNCOIndex(kSV)+(ii - 1)/fs/Tchip), CodeLen))).*...
%                     exp(2*pi*fd(kSV)*(1:1:NsamplesData)*dt*1i).*exp(phi(kSV)*1i);
%         end
%         %% filter and normalize SV signals before storing
% %         tic
%         x_delay2(kSV,:)  = filtfilt(h,1,x_delay1(kSV,:));
%         x_delay(kSV,:) = x_delay2(kSV,:)*sqrt((NsamplesData/sum(abs(x_delay2(kSV,:)).^2)));
% %         x1 = gpuArray(x_delay1(kSV,:));
% %         x2 = filtfilt(h,1,x1);
% %         x_delay(kSV,:) = gather(x2*sqrt((NsamplesData/sum(abs(x2).^2))));
% %         toc
%         %% filter and normalize s(t) for ZZB computation
%         s(kSV,:) = Code((1 + mod(round(PrevNCOIndex(kSV)+(ii - 1)/fs/Tchip), CodeLen)));
%         s(kSV,:) = filtfilt(h,1,s(kSV,:));
%         s(kSV,:) = s(kSV,:)*sqrt((NsamplesData/sum(s(kSV,:).^2)));
    %% generate C/A code, Doppler and carrier phase
    x_delay1(kSV,:) = circshift(caCode_filt(kSV,:), -round(PrevNCOIndex(kSV)*fs*Tchip));
    %% filter and normalize SV signals before storing
    x_delay2(kSV,:)  = x_delay1(kSV,:).*exp(2*pi*fd(kSV)*(1:1:NsamplesData)*dt*1i).*exp(phi(kSV)*1i);
    x_delay(kSV,:) = x_delay2(kSV,:)*sqrt((NsamplesData/sum(abs(x_delay2(kSV,:)).^2)));
    %% filter and normalize s(t) for ZZB computation
    s(kSV,:) = x_delay1(kSV,:);
    s(kSV,:) = s(kSV,:)*sqrt((NsamplesData/sum(s(kSV,:).^2)));
    end
%% gather outputs in a struct
% sigen.s = s;
sigen.x_delay = x_delay;
sigen.NsamplesLocal = NsamplesLocal;
sigen.NsamplesData = NsamplesData;





