function cn0 = estimateCn0(r,config,meanNoise)

%% load configuration file
eval(config)

%% find maximum value and its argument
[energy, maxPos] = max(r,[],2);


%% Estimate mean noise
if estimateTrueNoise == 0
    [~, r_samples] = size(r);  % Number of samples per signal
    chipSamples = ceil(Tc*fs); % Number of samples per chip


    low_range=maxPos-chipSamples;
    high_range=maxPos+chipSamples;
    % Remove autocorrelation from each signal (± 1 chip)
    r_clean=zeros(numSV,r_samples-chipSamples*2-1);
    for idx=1:numSV
        if low_range(idx)>1 && high_range(idx)<r_samples
            range = [1:low_range(idx)-1 high_range(idx)+1:r_samples];
        elseif low_range(idx)<=1
            range = high_range(idx)+1:(mod(low_range(idx)-2,r_samples)+1);
        else
            range = mod(high_range(idx),r_samples)+1:low_range(idx)-1;
        end
        r_clean(idx,:) = r(idx,range);
    end
    % Estimate mean noise
    meanNoise = mean(r_clean,2);
end

%% Estime snr and cn0
snr =(energy+meanNoise) ./ meanNoise;
cn0 = 10*log10(snr/(CodePeriod*CoherentIntegrations));


