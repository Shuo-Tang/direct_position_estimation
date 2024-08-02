function slope = generateCAFs_20ms_slope(usrInfo, satInfo, rawSignal, caCode, fd, corr)
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
satVel = satInfo(4:6, :)';
satCb = satInfo(7,:);
% satCd = satInfo(8,:);
% satEl = satInfo(9,:);
satPRN = satInfo(8, :);

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
% codeDiffms = [1856 7716 6377 18914] ./ 25375; 
%% memory allocation
r = 0;
%% DPE method
numSV = size(satInfo, 2);
nSamplesData = 25375;
%% generate local replica
rangeC = zeros(numSV, 1);
for ims = 1: 20
    rawSignal_1ms = rawSignal(:, 25375*(ims-1) + 1: 25375*ims);
    for kSV = 1%: numSV
        %% delay
        % corrcetion
        % caCodePRN = caCode(satPRN(kSV),:);
        code = generateCAcode(satPRN(kSV));
        geoRangeC = norm(satPos(kSV,:) - usrPos);
        % travelTime = geoRangeC/c;
        % corrSatPos = e_r_corr(travelTime, satPos(kSV,:)')';
        
        % compute delay
        % geoRangeC = norm(corrSatPos - usrPos);
        % iono = - 40.308 * 10e16 / (f0)^2 * ...
        %                 (1 - (6371/(6371 + 450)*sind(90 - satEl(kSV)))^2)^(-1/2);
        % trop = tropo(sin(satEl(kSV) * pi/180), 0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
    
        rangeC(kSV) = geoRangeC + c*(usrCb - satCb(kSV)) - corr(kSV, 1);%+ c*codeDiffms(kSV);% + corr(kSV, 1);% + trop + iono; %- 68.802*1e-3*c ; 
        fracDelayC = mod(rangeC(kSV)/c,CodePeriod); % - travelTimeOffset
        prevNCOIndexC = -fracDelayC/Tc;
        % prevNCOIndexC = 0;
        ii = 1: nSamplesData;
        x_shifted = code((1 + mod(round(prevNCOIndexC+(ii - 1)/fs/Tchip), length(code))));
        % IndDelay = fracDelayC*fs;
        % x_shifted = [caCode_10ms(ceil(fs*1e-3-IndDelay):end-1),...
        %     caCode_10ms( 1:ceil(fs*1e-3-IndDelay))];
        %% Doppler
        u = (usrPos-satPos(kSV,:))./norm(usrPos-satPos(kSV,:));
        fdC = -(usrVel-satVel(kSV,:))*u'*f0/c * (1+usrCd) + IF;
        % fdC = 421940.041529097;
        % fdC = fd(kSV);
        % fdC = fd(kSV);
        time = (0:25375) ./ 25375e3;
        DopplerMplC = exp((-2*pi*fdC*time(1:25375))*1i); % + freqOffset
        %% phase
        rangeCarrier =  geoRangeC + c*(usrCb - satCb(kSV));% + trop - iono;
        phiC = 2*pi*(rangeCarrier/lambda - corr(kSV, 2));
        phiMplC = exp(phiC*1i);
        %% combined 3 components
        % DopplerMplC = 1;
        % phiC = 0;
        % x_local = x_shifted.*DopplerMplC.*exp(phiC*1i);
        x_local = x_shifted;
        %% perform correlation
        r = r + real(sum(rawSignal_1ms(kSV,:).*DopplerMplC.*phiMplC.*conj(x_local)))^2;
        % z = abs(ifft(fft(DopplerMplC.*rawSignal_1ms(kSV,:)).*conj(fft(x_local)))).^2;
        caf_1ms(ims, kSV) = sum(rawSignal_1ms(kSV,:).*DopplerMplC.*phiMplC.*conj(x_local));
        % plot(z)
        % hold on
    end
end

%% test phase
phase_wrap = angle(caf_1ms);
phase = zeros(size(phase_wrap));
phase(1,:) = phase_wrap(1,:);
for kSV = 1%: numSV
    phase_add = 0;
    for ims = 2:20
        if abs(phase_wrap(ims-1, kSV) - phase_wrap(ims, kSV)) > pi
            phase_add = phase_add + 2*pi*sign(phase_wrap(ims-1, kSV) - phase_wrap(ims, kSV));
        end
        phase(ims, kSV) = phase_wrap(ims, kSV) + phase_add;
    end
end
sv_plot = 1;
figure, plot(phase_wrap(:,sv_plot))
hold on
plot(phase(:,sv_plot))
legend("wrapped phase", "phase")
slope = (phase(20,:) - phase(1,:))' ./ 19 * 1e3;
1;
% figure(101);
% title(sprintf("jump = %d ms", jumpMs))
% saveas(h1, sprintf("figs/jumpMsTest/jump_%dms", jumpMs))
% close(h1)





