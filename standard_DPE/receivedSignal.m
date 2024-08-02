function x = receivedSignal(sigen,config,CNo)


%% load configuration file
eval(config)

x_delay = sigen.x_delay;
NsamplesData = sigen.NsamplesData;


%% memory allocation
x_delay_noise=zeros(numSV,NsamplesData);



%% Add AWGN noise to the transmitted signals
noise = ( sqrt(1/2)*randn(1,NsamplesData) +1i* sqrt(1/2)*randn(1,NsamplesData) );
for kSV=1:numSV
    if CNo(kSV)<100
        % Sets amplitude assuming complex-noise power equal to 1
        %For CNo >=100 no noise is added.
        
        A       = sqrt(10^(CNo(kSV)/10)/fs);
        x_delay_noise(kSV,:) = A * x_delay(kSV,:);
    else
        x_delay_noise(kSV,:)=x_delay(kSV,:);
    end
end
%% Add noice to received signal and filter
x = sum(x_delay_noise, 1)+noise;

wn=pi*fn/(pi*fs/2);
h=fir1(order,wn);
x  = filtfilt(h,1,x);
