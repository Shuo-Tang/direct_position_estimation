function r = generateCAFs(config, rawSignal, usrInfo, satInfo, caCode)
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

%% number of samples calculation
NsamplesLocal=CodePeriod*fs*CoherentIntegrations;   % Number of samples of the Local Replica
NsamplesData=CodePeriod*fs*CoherentIntegrations*NonCoherentIntegrations;    %Number of samples of the Received Signal (Data).

%% generate local replica
r = 0;
for iSat = 1: nSat
    Code = caCode(satPRN(iSat),:);
    %% delay
    GeoRange = norm(satPos(iSat,:) - usrPos);
    Range = GeoRange + c*(usrCb - satCb(iSat));
    FracDelay = mod(Range/c, CodePeriod);
    PrevNCOIndex = -FracDelay/Tc;
    ii = 1: NsamplesData;
    x_shifted = Code(1 + mod(round(PrevNCOIndex+(ii - 1)/fs/Tchip), length(Code)));
    %% Doppler
    UnitDirection = (usrPos-satPos(iSat,:))./...
        vecnorm(usrPos-satPos(iSat,:),2,2);
    fd = -(usrVel-satVel(iSat,:))*UnitDirection'*f0/c...
        *(1+usrCd);
    Doppler = exp(2*pi*fd*(1:1:NsamplesLocal)*dt*1i);
    %% combine
    x_local= x_shifted.*Doppler;
    %% Perform Correlation
    r = r + abs(sum(rawSignal.*conj(x_local)))^2;
    % z = real(ifft(fft(rawSignal).*conj(fft(x_local)))).^2;
    % plot(z)
    % hold on
end
% 1;
% figure(101);
% title(sprintf("jump = %d ms", jumpMs))
% saveas(h1, sprintf("figs/jumpMsTest/jump_%dms", jumpMs))
% close(h1)





