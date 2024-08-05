function meanNoise = computeMeanNoise(config,sigen)

%% load configuration file
eval(config)

NsamplesData = sigen.NsamplesData;
if estimateTrueNoise == 0
    meanNoise = nan;
else
    %% Add AWGN noise to the transmitted signals
    mean_noise=zeros(1,Nexpe);
    for exp_idx=1:Nexpe
        noise = ( sqrt(1/2)*randn(1,NsamplesData) +1i* sqrt(1/2)*randn(1,NsamplesData) );
        r = correlateSignal(sigen,config,noise);
        mean_noise(exp_idx) = mean(mean(r,2));
    end
    meanNoise=mean(mean_noise);
end
end