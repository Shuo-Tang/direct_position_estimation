function [channel] = generateChannelInfo(prn, phase, delay, settings)
%GENERATECHANNELINFO generate channel for tracking 
%% Initialize all channels ================================================
channel                 = [];   % Clear, create the structure

channel.PRN             = 0;    % PRN number of the tracked satellite
channel.acquiredFreq    = 0;    % Used as the center frequency of the NCO
channel.codePhase       = 0;    % Position of the C/A  start

channel.status          = '-';  % Mode/status of the tracking channel
                                % "-" - "off" - no signal to track
                                % "T" - Tracking state
fs = settings.samplingFreq;
samplesPerCode = round(settings.samplingFreq / ...
                        (settings.codeFreqBasis / settings.codeLength));
%--- Copy initial data to all channels ------------------------------------
channel = repmat(channel, 1, settings.numberOfChannels);

%% Load information about each satellite --------------------------------
for ii = 1:settings.numberOfChannels
    channel(ii).PRN          = prn(ii);
    channel(ii).acquiredFreq = phase(ii);
    channel(ii).codePhase    = rem(round(delay(ii) * fs), samplesPerCode) + 1;
    % Set tracking into mode (there can be more modes if needed e.g. pull-in)
    channel(ii).status       = 'T';
end


end

