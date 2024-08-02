function acqResults = acquisition_10ms(longSignal, settings)
%Function performs cold start acquisition on the collected "data". It
%searches for GPS signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = acquisition(longSignal, settings)
%
%   Inputs:
%       longSignal    - 11 ms of raw signal from the front-end 
%       settings      - Receiver settings. Provides information about
%                       sampling and intermediate frequencies and other
%                       parameters including the list of the satellites to
%                       be acquired.
%   Outputs:
%       acqResults    - Function saves code phases and frequencies of the 
%                       detected signals in the "acqResults" structure. The
%                       field "carrFreq" is set to 0 if the signal is not
%                       detected for the given PRN number. 
 
%% Initialization =========================================================

% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
                        (settings.codeFreqBasis / settings.codeLength));

signal0DC1 = longSignal(1:(samplesPerCode*11)) - mean(longSignal(1:(samplesPerCode*11)));   %%Problems here....
signal0DC2 = longSignal((samplesPerCode*11 +1):(samplesPerCode*22)) - mean(longSignal((samplesPerCode*11+1):(samplesPerCode*22)));   %%Problems here....

% Find sampling period
ts = 1 / settings.samplingFreq;

% Find phase points of the local carrier wave 
phasePoints = (0 : (samplesPerCode-1)) * 2 * pi * ts;

% Number of the frequency bins for the given acquisition band (500Hz steps)
numberOfFrqBins = round(settings.acqSearchBand * 2) + 1;

% Generate all C/A codes and sample them according to the sampling freq.
caCodesTable = makeCaTable(settings); 


%--- Initialize arrays to speed up the code -------------------------------
% Search results of all frequency bins and code shifts (for one satellite)
results     = zeros(numberOfFrqBins, samplesPerCode);

% Carrier frequencies of the frequency bins
frqBins     = zeros(1, numberOfFrqBins);


%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 32);
% C/A code phases of detected signals
acqResults.codePhase    = zeros(1, 32);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 32);

fprintf('(');

% Perform search for all listed PRN numbers ...
for PRN = settings.acqSatelliteList
    signal_sample_index=0;
    %--- Perform DFT of C/A code ------------------------------------------
    caCodeFreqDom = conj(fft(caCodesTable(PRN, :)));
    for DWELL=1:settings.non_coherent_ms
        %% read one PRN block of signal
        input_signal_IF = longSignal(signal_sample_index+1 : (signal_sample_index+samplesPerCode));
        signal_sample_index=signal_sample_index+samplesPerCode;
        %% Correlate signals ======================================================
        %--- Make the correlation for whole frequency band (for all freq. bins)
        results=zeros(numberOfFrqBins,samplesPerCode);
        for frqBinIndex = 1:numberOfFrqBins
            %--- Generate carrier wave frequency grid (0.5kHz step) -----------
            frqBins(frqBinIndex) = settings.IF - ...
                (settings.acqSearchBand/2) * 1000 + ...
                0.5e3 * (frqBinIndex - 1);
            %--- Adjust the PRN delay drift due to code Doppler (see TSUI p.247)
%             delta_t=(DWELL-1)*1e-3;
%             shift_t_d=(frqBins(frqBinIndex)/settings.GPS_L1)*delta_t;
%             uniqFftPts = ceil((samplesPerCode) / 2);
%             fftFreqBins = (0 : uniqFftPts-1) * settings.samplingFreq/samplesPerCode;
%             fftFreqBins=[fftFreqBins -fftFreqBins(end:-1:1)];
%             
%             shift_fft=exp(1i*2*pi*fftFreqBins*shift_t_d);
%             caCodeFreqDom_shifted=caCodeFreqDom.*shift_fft;
            
            %--- Generate local sine and cosine -------------------------------
            sigCarr = exp(-1i*frqBins(frqBinIndex) * phasePoints);            
            %--- "Remove carrier" from the signal -----------------------------
            input_signal_BB=input_signal_IF.*sigCarr;            
            %--- Convert the baseband signal to frequency domain --------------
            input_signal_BB_fft = fft(input_signal_BB);           
            %--- Multiplication in the frequency domain (correlation in time
            %domain)
            %convCodeIQ = input_signal_BB_fft .* caCodeFreqDom_shifted;
            convCodeIQ = input_signal_BB_fft .* caCodeFreqDom;
            %--- Perform inverse DFT and store correlation results ------------
            acqRes = abs(ifft(convCodeIQ)) .^ 2;
            results(frqBinIndex, :) = results(frqBinIndex, :)+acqRes;
        end % frqBinIndex = 1:numberOfFrqBins
    end % frqBinIndex = 1:numberOfFrqBins
    
    %% plot the grid
    %Td=0:1:(samplesPerCode-1);
    %figure;
    %surf(Td,frqBins,results);

%% Look for correlation peaks in the results ==============================
    % Find the highest peak and compare it to the second highest peak
    % The second peak is chosen not closer than 1 chip to the highest peak
    
    %--- Find the correlation peak and the carrier frequency --------------
    [peakSize frequencyBinIndex] = max(max(results, [], 2));
    
    % --- plot the time delay line at the maximum peak
    %figure;
    %stem(results(frequencyBinIndex,:));

    %--- Find code phase of the same correlation peak ---------------------
    [peakSize codePhase] = max(max(results));

    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
    samplesPerCodeChip   = round(settings.samplingFreq / settings.codeFreqBasis);
    excludeRangeIndex1 = codePhase - samplesPerCodeChip;
    excludeRangeIndex2 = codePhase + samplesPerCodeChip;

    %--- Correct C/A code phase exclude range if the range includes array
    %boundaries
    if excludeRangeIndex1 < 2
        codePhaseRange = excludeRangeIndex2 : ...
                         (samplesPerCode + excludeRangeIndex1);
                         
    elseif excludeRangeIndex2 >= samplesPerCode
        codePhaseRange = (excludeRangeIndex2 - samplesPerCode) : ...
                         excludeRangeIndex1;
    else
        codePhaseRange = [1:excludeRangeIndex1, ...
                          excludeRangeIndex2 : samplesPerCode];
    end

    %--- Find the second highest correlation peak in the same freq. bin ---
    secondPeakSize = max(results(frequencyBinIndex, codePhaseRange));

    %--- Store result -----------------------------------------------------
    acqResults.peakMetric(PRN) = peakSize/secondPeakSize;
    
    % If the result is above threshold, then there is a signal ...
    if (peakSize/secondPeakSize) > settings.acqThreshold

%% Fine resolution frequency search =======================================
        
        %--- Indicate PRN number of the detected signal -------------------
        fprintf('%02d ', PRN);
        
        %--- Generate 10msec long C/A codes sequence for given PRN --------
        caCode = generateCAcode(PRN);
        
        codeValueIndex = floor((ts * (1:10*samplesPerCode)) / ...
                               (1/settings.codeFreqBasis));
                           
        longCaCode = caCode((rem(codeValueIndex, 1023) + 1));
    
        %--- Remove C/A code modulation from the original signal ----------
        % (Using detected C/A code phase)
        xCarrier1 = ...
            signal0DC1(codePhase:(codePhase + 10*samplesPerCode-1)) ...
            .* longCaCode;
        
        xCarrier2 = ...
            signal0DC2(codePhase:(codePhase + 10*samplesPerCode-1)) ...
            .* longCaCode;
        %--- Compute the magnitude of the FFT, find maximum and the
        %associated carrier frequency
        
        %--- Find the next highest power of two and increase by 8x --------
        fftNumPts = 8*(2^(nextpow2(length(xCarrier1))));
        
        %--- Compute the magnitude of the FFT, find maximum and the
        %associated carrier frequency 
        fftxc1 = abs(fftshift(fft(xCarrier1, fftNumPts))); 
        fftxc2 = abs(fftshift(fft(xCarrier2, fftNumPts))); 
        
        [fftMax1, fftMaxIndex1] = max(fftxc1);
        [fftMax2, fftMaxIndex2] = max(fftxc2);
        
        fftFreqBins=settings.samplingFreq/2*(-1:2/fftNumPts:1-2/fftNumPts);
        
        %--- Compute the magnitude of the FFT, find maximum and the
        %associated carrier frequency 
        if (fftMax1>fftMax2)
            acqResults.carrFreq(PRN)  = fftFreqBins(fftMaxIndex1);
            if (settings.enable_acquisition_debug_plots==1)
                figure;
                plot(fftFreqBins,fftxc1);
            end
        else
            acqResults.carrFreq(PRN)  = fftFreqBins(fftMaxIndex2);
            if (settings.enable_acquisition_debug_plots==1)
                figure;
                plot(fftFreqBins,fftxc2);
            end
        end
         acqResults.carrFreq(PRN) =acqResults.carrFreq(PRN);
         acqResults.codePhase(PRN) = codePhase;
    
    else
        %--- No signal with this PRN --------------------------------------
        acqResults.carrFreq(PRN)=0;
        acqResults.codePhase(PRN)=0;
        fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold
    
end    % for PRN = satelliteList

%=== Acquisition is over ==================================================
fprintf(')\n');
