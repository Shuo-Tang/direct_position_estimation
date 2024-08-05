function [sigen] = signalGen(config)


%% load configuration file
eval(config)


%% number of samples calculation
NsamplesLocal=CodePeriod*fs*CoherentIntegrations;   % Number of samples of the Local Replica
NsamplesData=CodePeriod*fs*CoherentIntegrations*NonCoherentIntegrations;    %Number of samples of the Received Signal (Data).


%% memory allocation
Range=zeros(1,numSV);
x_local=zeros(numSV,NsamplesLocal);
fft_local=zeros(numSV,NsamplesLocal);
x_delay=zeros(numSV,NsamplesData);


%% Compute range and fractional delays for each SV

for kSV=1:numSV 
    Range(kSV)                               =   norm(SatPosition(kSV,:) - UserPosition);    
end

FracDelay=mod(Range/c,CodePeriod);


%% Generate local replica and delayed signals according to the computed delays
   
PrevNCOIndex    =  -  FracDelay/Tc;
randomDelay= 0;
for kSV=1:numSV
%     Code                                    =   genCAcode(SatPRN(kSV));
    Code = CAgenerator_L5IQ(SatPRN(kSV),'L5I');
    Tchip                                   =   CodePeriod / length(Code);
    ii                                      =   1 : NsamplesLocal;          % generating LGenBlocks samples
    x_local(kSV,:)                                 =   Code((1 + mod(round((ii - 1) / fs / Tchip), length(Code))));
%     fft_local(kSV,:) = fft(x_local(kSV,:),Nsamples);
    ii                                      =   1 : NsamplesData;
    x_delay(kSV,:)                                 =   Code((1 + mod(round(PrevNCOIndex(kSV)+randomDelay+(ii - 1) / fs / Tchip), length(Code))));
end

%% Filter local signal and generate its FFT 

wn=pi*fn/(pi*fs/2);
h=fir1(order,wn);
for kSV=1:numSV
    x_delay(kSV,:)  = filtfilt(h,1,x_delay(kSV,:));
    x_delay(kSV,:)  = x_delay(kSV,:)*sqrt(length(x_delay(kSV,:))/sum(abs(x_delay(kSV,:)).^2));
    x_local(kSV,:)  = filtfilt(h,1,x_local(kSV,:));
    x_local(kSV,:)  = x_local(kSV,:)*sqrt(length(x_local(kSV,:))/sum(abs(x_local(kSV,:)).^2));
    fft_local(kSV,:) = fft(x_local(kSV,:),NsamplesLocal);
end

%% Normalize Received Signal Power after filtering

% for kSV=1:numSV
%     x_delay(kSV,:)  = x_delay(kSV,:)*sqrt((NsamplesData/sum(x_delay(kSV,:).^2)));
% end


%% gather outputs in a struct
sigen.x_delay = x_delay;
sigen.x_local = x_local;
sigen.fft_local = fft_local;
sigen.randomDelay = randomDelay;
sigen.NsamplesLocal = NsamplesLocal;
sigen.NsamplesData = NsamplesData;



