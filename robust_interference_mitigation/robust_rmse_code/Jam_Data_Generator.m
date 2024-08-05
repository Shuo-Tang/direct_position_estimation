
function [y] = Jam_Data_Generator(config)
eval(config)
%% number of samples calculation
NsamplesData=CodePeriod*fs*CoherentIntegrations*NonCoherentIntegrations;    %Number of samples of the Received Signal (Data).

fOfjam=1e4;%1e3;
A = 1;
t = (1:1:NsamplesData)*dt;
jam1=A*cos(2*pi*t*fOfjam);
jam2=A*sin(2*pi*t*fOfjam);
y = jam1 + 1j.*jam2;
y = repmat(y, [numSV, 1]);
end
