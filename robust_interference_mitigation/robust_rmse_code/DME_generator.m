function [y_new] = DME_generator(config)
%% number of samples calculation

eval(config)
NsamplesData=CodePeriod*fs*CoherentIntegrations*NonCoherentIntegrations;    %Number of samples of the Received Signal (Data).
alpha = 4.5e11;
tao = 6e-6;
A = 1;
lambda = 2700;
f_c = 5e6;
y = zeros(1, NsamplesData);
num_time = NsamplesData/fs;
window_half = 800;
interval = 0;
t=(1:1:NsamplesData)*dt;
for i = 1 : num_time*lambda
    interval = interval+exprnd(1/lambda);
    per_ceil = ceil(interval * fs - window_half);
    per_floor = ceil(interval * fs + window_half);
    if(per_ceil <= 0 || per_floor > NsamplesData)
        1;
    else
    DME_per = exp(-alpha/2*(t(per_ceil:per_floor)-tao-interval).^2)...
        + exp(-alpha/2*(t(per_ceil:per_floor)+tao-interval).^2);
    y(per_ceil : per_floor)=y(per_ceil:per_floor) + DME_per;
    end
end
S_AM=A*exp(1j*2*pi*f_c*t);
y_new=y.*S_AM;
y_new = repmat(y_new, [numSV, 1]);
% DME_Carrier=[real(y_new);imag(y_new)];
end
