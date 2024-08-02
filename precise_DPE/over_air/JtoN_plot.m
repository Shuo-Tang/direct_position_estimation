addpath include             % The software receiver functions
addpath geoFunctions        % Position calculation related functions 
clear all
% Initialize parameters
settings = initSettings();
numOfA=2500;
i=100;
num=zeros(1,numOfA/i+1);
pow=zeros(1,numOfA/i+1);
CtoN=zeros(1,numOfA/i+1);
varD=zeros(1,numOfA/i+1);
varP=zeros(1,numOfA/i+1);

h=1; % Number of loops
for a=0:i:numOfA
    filename='/agilent_cap2.dat';
    fid = fopen(filename,'r');
    [JtoN0,rawdata_process]=jamming_data_produce(fid,a); % Get JtoN0
    nfid=fopen('newrawdata.dat','r');
    fseek(nfid, 8*settings.skipNumberOfBytes, 'bof'); 
    settings.samplingFreq       = 4e6;
    settings.codeFreqBasis      = 1.023e6;      
    settings.codeLength         = 1023;
    samplesPerCode = round(settings.samplingFreq/(settings.codeFreqBasis / settings.codeLength));
    [newdata, count] = fread(nfid, [2, 1000*samplesPerCode], settings.dataType);
    newdata = newdata(1,:) + newdata(2,:).*1i;
    %% Notch filter
%     newdata1=newdata;
%     nfilter = struct( 'z0', 0, ...              % filter zero
%                   'xi', 0, ...              % output of the AR part  
%                   'ka', 0.8, ...            % contraction factor
%                   'mu', 0 );
%               numEpochs=1;
% for ii = 1:numEpochs,
% 
%     if rem(ii, 100) == 0,
%         ii          % display some info on the processing status
%     end
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                               Read the data                             %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%    
%             
%         rawdata = newdata1;
%     
% %     rawdata = rawdata.';
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %                     Perform notch filtering                         %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [ filtdata, nfilter, z0 ] = adaptivenotch( rawdata, nfilter );
%     
%     
%     % Use the filtered samples only if the interference is clearly detected
%     if mean( abs( z0 ) ) > 0.75,
%        newdata=filtdata;
%     else
%        newdata=rawdata;
%     end
% end
    
    
    %% Robustness pre-correlation
%     fftNdata=fft(newdata);
%     RfftNdata=fftNdata./abs(fftNdata);
%     newdata=ifft(RfftNdata);
%     tfid=fopen('newrawdatafortracking.dat','wb');
%     Tnewdata=[real(newdata);imag(newdata)];
%     fwrite(tfid,Tnewdata,'float');
%     fclose(tfid);
%     tfid=fopen('newrawdatafortracking.dat','rb');

% Get acqResults
    acqResults = acquisition(newdata, settings);
    
% Get channel
     if (any(acqResults.carrFreq))
        channel = preRun(acqResults, settings);
        % Get trackResults
        [trackResults, channel,varianceDLL,variancePLL]= tracking(fid, channel, settings,a);

    else
        % No satellites to track, exit
%         trackResults = [];
         varianceDLL=0;
         variancePLL=0;
%         return;
     end
   
    
    fclose(fid);
% Get trackResults
%     [trackResults, channel]= tracking(fid, channel, settings);
%   fclose(tfid);
    varD(h)=varianceDLL;
    varP(h)=variancePLL;
    pow(h)=JtoN0;
    num(h)=length(nonzeros(acqResults.carrFreq));
    CtoN(h)=acqResults.CtoN(11);
    h=h+1;
    fclose(nfid);
end

% Get figure: SVs with J/N0 
figure;
plot(pow,num);
set(gca,'YLim',[0,9]);
title('Cosine Jamming signal for number of availabe SVs');
%title('Random Jamming signal toward Normalized');
xlabel('J/N0 (dB)');
ylabel('Number of available SVs');

% Get figure: C/N of satelattie 11 with J/N0
figure;
plot(pow,CtoN);
axis tight; 
title('Cosine Jamming signal for SNR');
xlabel('J/N0 (dB)');
ylabel('C/N for SV 11(dB)');

%Get figure: Discriminator with J/NO
figure;
plot(pow,varD,'--');
hold on
plot(pow,varP,'-');
% axis tight;
title('Cosine Jamming signal for Discriminator');
legend('VarianceDLL','VariancePLL');
xlabel('J/N0 (dB)');
ylabel('Discriminators for SV 11');

