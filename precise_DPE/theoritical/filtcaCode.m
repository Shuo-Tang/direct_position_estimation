function caCode_filt = filtcaCode(config)

eval(config)
%% number of samples calculation
NsamplesLocal=CodePeriod*fs*CoherentIntegrations;   % Number of samples of the Local Replica
NsamplesData=CodePeriod*fs*CoherentIntegrations*NonCoherentIntegrations;    %Number of samples of the Received Signal (Data).
%% setting for filter
fn=2e6;
order=36;
wn=pi*fn/(pi*fs/2);
h=fir1(order,wn);
%% memory set
caCode = zeros(numSV,NsamplesData);
caCode_filt = zeros(numSV,NsamplesData);
%% generate CAcode
for kSV = 1:numSV
    Code = genCAcode(SatPRN(kSV));
    ii = 1 : NsamplesData;
    %% generate C/A code, Doppler and carrier phase
    caCode(kSV,:) = Code((1 + mod(round((ii - 1)/fs/Tchip), CodeLen)));
    %% filter and normalize SV signals before storing
    caCode_filt(kSV,:)  = filtfilt(h,1,caCode(kSV,:));
end

