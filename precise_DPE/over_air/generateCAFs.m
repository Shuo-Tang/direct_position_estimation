function [r] = generateCAFs(usrInfo, satInfo, rawSignal, caCode, timeDiff)
% This function is designed for computing CAFs
% input:  usrInfo   -- receiver position, velocity, clock bias and clock
%                   drift information
%         satInfo   -- satellite position, velocity, clock bias, clock
%                   drift, el, az and PRN information
%         rawsignal -- received raw signal
%         caCode    -- C/A code generated for all GPS satellites based on
%                   sampling frequency
%         timeDiff  -- sampling moment difference w.r.t the minimum one 
%                      (if clock bias is searched, this will not be used)
%% load configuration file
config = 'ConfigFile';
eval(config)

satPos = satInfo(1:3, :)';
satVel = satInfo(5:7, :)';
satCb = satInfo(4,:);
% satCd = satInfo(8,:);
satEl = satInfo(9,:);
satPRN = satInfo(11, :);

usrPos = usrInfo(1:3);
usrVel = usrInfo(5:7);
usrCb = usrInfo(4);
usrCd = usrInfo(8);
IF = usrInfo(9);
% freqOffset = usrInfo.freqOffset;

% nExtraCode = length(rawSignal) - NsamplesLocal;
% if nExtraCode < 0
%     rawSignal = [rawSignal, zeros(1, -nExtraCode)];
% end

%% memory allocation
r = 0;
%% DPE method
numSV = size(satInfo, 2);
nSamplesData = length(rawSignal);
%% generate local replica
for kSV = 1:numSV
    %% delay
    % corrcetion
    geoRangeC = norm(satPos(kSV,:) - usrPos);
    travelTime = geoRangeC/c;
    corrSatPos = e_r_corr(travelTime, satPos(kSV,:)')';
    trop = tropo(sin(satEl(kSV) * pi/180), 0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
    iono = 40.308 * 34e16 / (1575.42e6)^2; 
    % compute delay
    geoRangeC = norm(corrSatPos - usrPos);
    rangeC = geoRangeC + c*(usrCb - satCb(kSV)) + trop - timeDiff(kSV)*c - iono;% - 68.802*1e-3*c ; 
    fracDelayC = mod(rangeC/c,CodePeriod); % - travelTimeOffset
    % prevNCOIndexC = -fracDelayC/Tc;
    % ii = 1: nSamplesData;
    % x_shifted = code((1 + mod(round(prevNCOIndexC+(ii - 1)/fs/Tchip), length(code))));
    IndDelay = fracDelayC*fs;
    caCodePRN = caCode(satPRN(kSV),:);
    x_shifted = [caCodePRN(ceil(fs*1e-3-IndDelay):end-1),...
        caCodePRN( 1:ceil(fs*1e-3-IndDelay))];
    %% Doppler
    u = (usrPos-corrSatPos)./vecnorm(usrPos-corrSatPos,2,2);
    fdC = -(usrVel-satVel(kSV,:))*u'*f0/c * (1+usrCd) + IF;
    DopplerMplC = exp(2*pi*fdC*(1:1:nSamplesData)*dt*1i); % + freqOffset
    %% phase
    phiC = 2*pi*((rangeC)/lambda);
    %% combined 3 components
    % DopplerMplC = 1;
    % phiC = 0;
    x_local = x_shifted.*DopplerMplC.*exp(phiC*1i);

    %% perform correlation
    r = r + abs(sum(rawSignal(kSV, :).*conj(x_local)))^2;
end





