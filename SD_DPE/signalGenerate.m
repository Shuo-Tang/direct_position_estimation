function [x_tx] = signalGenerate(config, usrInfo, satInfo, caCode)


%% load configuration file and data
eval(config)
satPos = satInfo(:, 1:3);
satCb = satInfo(:, 4);
satVel = satInfo(:, 5:7);
satCd = satInfo(:, 8);
satPRN = satInfo(:, 9);
nSat = length(satPRN);

usrPos = usrInfo(1:3);
usrCb = usrInfo(4);
usrVel = usrInfo(5:7);
usrCd = usrInfo(8);

% setting for filter
% fn=2e6;
% order=36;
% wn=pi*fn/(pi*fs/2);
% h=fir1(order,wn);
%% number of samples calculation
NsamplesLocal=CodePeriod*fs*CoherentIntegrations;   % Number of samples of the Local Replica
NsamplesData=CodePeriod*fs*CoherentIntegrations*NonCoherentIntegrations;    %Number of samples of the Received Signal (Data).

%% memory allocation
xTransmit = zeros(nSat, NsamplesData);

%% generate local replica
for kSV = 1:nSat
    Code = caCode(satPRN(kSV),:);
    %% delay
    GeoRange = norm(satPos(kSV,:) - usrPos);
    Range = GeoRange + c*(usrCb - satCb(kSV));
    FracDelay = mod(Range/c, CodePeriod);
    PrevNCOIndex = -FracDelay/Tc;
    ii = 1: NsamplesData;
    x_shifted = Code(1 + mod(round(PrevNCOIndex+(ii - 1)/fs/Tchip), length(Code)));
    %% Doppler
    UnitDirection = (usrPos-satPos(kSV,:))./...
        vecnorm(usrPos-satPos(kSV,:),2,2);
    fd = -(usrVel-satVel(kSV,:))*UnitDirection'*f0/c...
        *(1+usrCd);
    Doppler = exp(2*pi*fd*(1:1:NsamplesLocal)*dt*1i);
    %% combine
    xTransmit(kSV, :) = x_shifted.*Doppler;
    %% filter and normalize the transmission signal
    % xFilter = filtfilt(h,1,xLocal);
end
x_tx = sum(xTransmit, 1);
end



