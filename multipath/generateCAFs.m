function r = generateCAFs(usrInfo, satInfo, rawSignal, caCode)
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
%% Load Predefined parameters
c = 299792458;          % Speed of light.
CodePeriod = 1e-3;
fs = 20.46e6;

satPos = satInfo(:,1:3);
satVel = satInfo(:,4:6);
satCb = 0;
satPRN = satInfo(:, 9);

usrPos = usrInfo(1:3);
usrCb = usrInfo(4);
usrVel = usrInfo(5:7);
usrCd = usrInfo(8);

nSat = length(satPRN);

%% generate local replica
r = 0;
range = zeros(nSat, 1);
for iSat = 1: nSat
    %% Delay Shift
    geoRange = norm(satPos(iSat,:) - usrPos);
    range(iSat) = geoRange + c* usrCb;
    fracDelayC = mod(range(iSat) / c, CodePeriod);
    caCodePRN = caCode(iSat,:);    
    IndDelay = fracDelayC * fs;

    % ii = 1: nSamplesData;
    % prevNCOIndexC = -fracDelayC/Tc;
    % x_shifted = code((1 + mod(round(prevNCOIndexC+(ii - 1)/fs/Tchip), length(code))));
    % x_shifted1 = [caCodePRN(round(fs*1e-3-IndDelay):end-1),...
    %     caCodePRN( 1:round(fs*1e-3-IndDelay))];
    x_local = circshift(caCodePRN, round(IndDelay));
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





