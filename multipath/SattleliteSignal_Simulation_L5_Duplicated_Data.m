clear all
fs=40e6;
ts = 1/fs;   % Sampling period in sec
numOfinterval=3e3;
N=fs*1e-3*numOfinterval;
tc = 1/(1.023e6*10);  % C/A chip period in sec
tb=1/50;
codeLength = 1023*10;
samplesPerCode=fs*codeLength*tc;
f_c=1300;       % Doppler Shift
pfid=fopen('sate_11_L5_Duplicated.dat','wb');
loops=1;
Mflag = true;
% addpath('D:\MultiLayerPerceptron\LandMobileMultipathChannelModel32')
MultipathChannel_settings;
for iloop=1:loops
codeValueIndex = ceil((ts * ((iloop-1)*samplesPerCode*numOfinterval+1:...
    iloop*samplesPerCode*numOfinterval)) / tc);
codeValueIndex=rem(codeValueIndex,numOfinterval*codeLength);
codeValueIndex(codeValueIndex==0)=(numOfinterval*codeLength);
bitValueIndex = ceil((ts * ((iloop-1)*samplesPerCode*numOfinterval+1:...
    iloop*samplesPerCode*numOfinterval)) / tb);
bitValueIndex=rem(bitValueIndex,numOfinterval/20);
bitValueIndex(bitValueIndex==0)=numOfinterval/20;
t=((iloop-1)*N+1:iloop*N)/fs;
yecho = 0;

%% Bit Information Generator
Bit=sign(randn(1,numOfinterval/20));
Bit_sample=Bit(bitValueIndex);
%% C/A code Generator for Sate 11

[CA_I]=CAgenerator_L5IQ(11,'L5I');
[CA_Q]=CAgenerator_L5IQ(11,'L5Q');
%% White Noise Generator
% noise=randn(1,samplesPerCode*numOfinterval)+randn(1,samplesPerCode*numOfinterval)*1i;
% noise=sqrt(10^(SNR/10)).*noise;

%% Signal combination
Carrier=exp(1j*2*pi*(f_c)*t);
CAcode_I=repmat(CA_I,1,numOfinterval);
CAcode_Q=repmat(CA_Q,1,numOfinterval);%+CAcode_Q(codeValueIndex).*1i
y=Bit_sample.*(CAcode_I(codeValueIndex)).*Carrier;%+noise;
z = y;
skip = 300*samplesPerCode;
if (iloop>1)
    skip = 0;
end
skip = 0;
Parameters.NumberOfSteps = (1e-3*numOfinterval-skip/fs)*Parameters.SampFreq;
if Mflag
    for dhv=1:Parameters.NumberOfSteps
        TimeVec(end)=dhv/Parameters.SampFreq;
        % --- "drunken" pedestrian movement example ---
        ActualSpeed=Parameters.MaximumSpeed/2/3.6*(1+sin(TimeVec(end)));   
        SpeedVec(dhv)=ActualSpeed;              % m/s
        ActualHeading=20*sin(TimeVec(end));     % Deg (North == 0, East == 90, South == 180, West == 270)
        % --- generate CIR ---
        [TheChannelObject,LOS,LOSDelays,ComplexOutputVec,DelayVec,EchoNumberVec,WayVec(dhv),TimeVec(dhv)]=...
        generate(TheChannelObject,ActualSpeed,ActualHeading,Parameters.SatElevation,Parameters.SatAzimut);
        if isempty(DelayVec)
            interval = 4e-7;
        else
            interval = max(max(DelayVec)+1e-7,4e-7);
        end
        NLOSdelay = (1:interval*fs)/fs;
        ind = round(DelayVec*fs);
        index  = unique(ind);
        NLOSloss = zeros(1,length(NLOSdelay));
        for i = 1:length(index)
            NLOSloss(index(i)+1) = sum(ComplexOutputVec(ind == index(i)));
        end
        LOSdelay = (1:interval*fs)/fs;
        LOSloss = zeros(1,length(LOSdelay));
        ind = round(LOSDelays*fs);
        index  = unique(ind);
        for i = 1:length(index)
            LOSloss(index(i)+1) = sum(LOS(ind == index(i)));
        end
        LOSloss(1) = 1;
        step = fs/Parameters.SampFreq;
        ytotal = cconv(y(((dhv-1)*step + 1:dhv*step)+skip),LOSloss) ...
            + [yecho, zeros(1,step+length(NLOSloss)-length(yecho)-1)];% ...
            % + cconv(y((dhv-1)*step + 1:dhv*step),NLOSloss);%
        yecho = ytotal(step+1:end);
        y(((dhv-1)*step + 1:dhv*step)+skip) = ytotal(1:step);
        
    end%
end
    
acum=1;
for size=1:length(codeValueIndex)/1e7
    sate1=y(acum:size*1e7);
    sate=[real(sate1);imag(sate1)];
    acum=acum+1e7;
    fwrite( pfid, sate, 'short' );
end

end
fclose(pfid);

