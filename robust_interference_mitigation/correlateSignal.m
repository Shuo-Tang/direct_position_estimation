function r = correlateSignal(sigen,config,x_delay_noise)


%% load configuration file
eval(config)

fft_local = sigen.fft_local;
NsamplesLocal = sigen.NsamplesLocal;


%% memory allocation
r=zeros(numSV,NsamplesLocal);
    %% Perform NonCoherentIntegrations times non coherent integrations of CoherentIntegrations times coherent integrations.
for kSV=1:numSV
    for idx_nc=1:NonCoherentIntegrations
        r(kSV,:)=r(kSV,:)+abs(ifft(fft(x_delay_noise(1,NsamplesLocal*(idx_nc-1)+1:NsamplesLocal*idx_nc),NsamplesLocal) .* conj(fft_local(kSV,:)))).^2;
    end
end
end

