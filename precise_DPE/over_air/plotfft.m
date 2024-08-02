clear all
filename='/agilent_cap2.dat';
fid = fopen(filename,'r');
fs=10e6;
fd=1.023e6; 
cl=1023;
fj=100;
samplesPerCode = round(fs/(fd / cl));

n=0:2:100*2*samplesPerCode-1;
jam1=1000*cos(2*pi*fj*n);
m=1:2:100*2*samplesPerCode-1;
jam2=1000*cos(2*pi*fj*m);
jam=[jam1;jam2];
[data, count] = fread(fid, [2, 100*samplesPerCode], 'float');
data=data+jam;
data = data(1,:) + data(2,:).*1i;

[sigspec,freqv]=pwelch(data, 32758, 2048, 16368, fs,'twosided');
figure;
plot(([-(freqv(length(freqv)/2:-1:1));freqv(1:length(freqv)/2)])/1e6, ...
            10*log10([sigspec(length(freqv)/2+1:end);sigspec(1:length(freqv)/2)]));
        
  
  