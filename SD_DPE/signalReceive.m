function xReceived = signalReceive(config, transmitSignal)


%% load configuration file
eval(config)

numSV = size(transmitSignal, 1);
NsamplesData = size(transmitSignal, 2);
CNo = eachCNo * ones(numSV, 1);
% set filter config
order = 36;
fn = 2e6;
wn = pi*fn/(pi*fs/2);
h = fir1(order,wn);
%% memory allocation
xNoised = zeros(numSV,NsamplesData);

%% Add AWGN noise to the transmitted signals
for kSV=1:numSV       
    A = sqrt(10^(CNo(kSV)/10)/fs);
    noise = sqrt(1/2)*randn(1,NsamplesData) + 1i*sqrt(1/2)*randn(1,NsamplesData);
%     noise = 0;
    xNoised(kSV,:) = A*transmitSignal(kSV,:) + noise;
    xNoised(kSV,:) = filtfilt(h,1,xNoised(kSV,:));
end
xReceived = sum(xNoised, 1);

