clear all

addpath('D:\GPSL5simulator\include')
addpath('D:\GPSL5simulator\geoFunctions')
settings = initSettings();
samplesPerCode = round(settings.samplingFreq / ...
                           (settings.codeFreqBasis / settings.codeLength));
CN0 = 40;
mode_linear = true;
MLPpath = 'D:\MultiLayerPerceptron\MLP_101_noise.mat';
SNRlin = 10^((CN0-10*log10(settings.samplingFreq))/10);
ts = 1/settings.samplingFreq;
tc = 1/settings.codeFreqBasis; 
codeValueIndex = ceil((ts * (1:samplesPerCode)) / tc);                                             
[fid, message] = fopen(settings.fileName, 'rb');   
dataAdaptCoeff = settings.fileType;
fseek(fid, dataAdaptCoeff*settings.skipNumberOfBytes, 'bof'); 
[data, count] = fread(fid, [2, 22*samplesPerCode], settings.dataType);
data = data(1,:) + data(2,:).*1i; %Inphase and Quadrature
PowS = sum(abs(data).^2)/length(data);
PowN = PowS/SNRlin; 
noise = sqrt(PowN/2).* (randn(1,length(data)) + 1j* randn(1,length(data)));
data = data + noise;
% acqResults = acquisition(data, settings);
% 
% if (any(acqResults.carrFreq))
%         channel = preRun(acqResults, settings);
%         showChannelStatus(channel, settings);
% else
%     % No satellites to track, exit
%     disp('No GNSS signals detected, signal processing finished.');
%     trackResults = [];
%     return;
% end
load('Channel_information.mat')
% channel.acquiredFreq = 1300;
[trackResults, channel] = tracking(fid, channel, settings, PowN, mode_linear, MLPpath);
plotTracking(1:settings.numberOfChannels, trackResults, settings);
fclose(fid);