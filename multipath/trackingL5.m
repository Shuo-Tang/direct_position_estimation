function [trackResults, channel]= tracking(fid, channel, settings, PowN, mode_linear, MLPpath)
% Performs code and carrier tracking for all channels.
%
%[trackResults, channel] = tracking(fid, channel, settings)
%
%   Inputs:
%       fid             - file identifier of the signal record.
%       channel         - PRN, carrier frequencies and code phases of all
%                       satellites to be tracked (prepared by preRum.m from
%                       acquisition results).
%       settings        - receiver settings.
%   Outputs:
%       trackResults    - tracking results (structure array). Contains
%                       in-phase prompt outputs and absolute spreading
%                       code's starting positions, together with other
%                       observation data from the tracking loops. All are
%                       saved every millisecond.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Based on code by DMAkos Oct-1999
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: tracking.m,v 1.14.2.31 2006/08/14 11:38:22 dpl Exp $

%% Initialize result structure ============================================

% Channel status
trackResults.status         = '-';      % No tracked signal, or lost lock

% The absolute sample in the record of the C/A code start:
trackResults.absoluteSample = zeros(1, settings.msToProcess);
trackResults.PRN_start_sample=zeros(1, settings.msToProcess);

% Freq of the C/A code:
trackResults.codeFreq       = inf(1, settings.msToProcess);

% Frequency of the tracked carrier wave:
trackResults.carrFreq       = inf(1, settings.msToProcess);

% Outputs from the correlators (In-phase):
trackResults.I_P            = zeros(1, settings.msToProcess);
trackResults.I_E            = zeros(1, settings.msToProcess);
trackResults.I_L            = zeros(1, settings.msToProcess);

% Outputs from the correlators (Quadrature-phase):
trackResults.Q_E            = zeros(1, settings.msToProcess);
trackResults.Q_P            = zeros(1, settings.msToProcess);
trackResults.Q_L            = zeros(1, settings.msToProcess);

% Loop discriminators
trackResults.dllDiscr       = inf(1, settings.msToProcess);
trackResults.dllDiscrFilt   = inf(1, settings.msToProcess);
trackResults.pllDiscr       = inf(1, settings.msToProcess);
trackResults.pllDiscrFilt   = inf(1, settings.msToProcess);

%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Initialize tracking variables ==========================================

PRN_start_sample = 0; %javi
codePeriods = settings.msToProcess;     % For GPS one C/A code is one ms

%--- DLL variables --------------------------------------------------------
% Define early-late offset (in chips)
earlyLateSpc = settings.dllCorrelatorSpacing;

% Summation interval
PDIcode = 0.001;

% Calculate filter coefficient values
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
                                    settings.dllDampingRatio, ...
                                    1.0);

%--- PLL variables --------------------------------------------------------
% Summation interval
PDIcarr = 0.001;

% Calculate filter coefficient values
[tau1carr, tau2carr] = calcLoopCoef(settings.pllNoiseBandwidth, ...
                                    settings.pllDampingRatio, ...
                                    0.25);
hwb = waitbar(0,'Tracking...');

if (settings.fileType==1) 
    dataAdaptCoeff=1;
end
if (settings.fileType==2) 
    dataAdaptCoeff=2;
end

if (settings.fileType==4) 
    dataAdaptCoeff=4;
end

if (settings.fileType==8) 
    dataAdaptCoeff=8;
end




%% Start processing channels ==============================================
    for channelNr = 1:settings.numberOfChannels
        % Only process if PRN is non zero (acquisition was successful)
        if (channel(channelNr).PRN ~= 0)
            % Save additional information - each channel's tracked PRN
            trackResults(channelNr).PRN     = channel(channelNr).PRN;

            % Move the starting point of processing. Can be used to start the
            % signal processing at any point in the data record (e.g. for long
            % records). In addition skip through that data file to start at the
            % appropriate sample (corresponding to code phase). Assumes sample
            % type is schar (or 1 byte per sample) 
              fseek(fid, ...
                dataAdaptCoeff*(settings.skipNumberOfBytes + channel(channelNr).codePhase-1), ...
                'bof');


            % Get a vector with the C/A code sampled 1x/chip
            [CA_I]=CAgenerator_L5IQ(channel(channelNr).PRN,'L5I');
            [CA_Q]=CAgenerator_L5IQ(channel(channelNr).PRN,'L5Q');
            caCode = conj(CA_I );%+ CA_Q.*1i);
            % Then make it possible to do early and late versions
            if mode_linear
                caCode = [caCode(end) caCode caCode(1)];
            else
                caCode = [caCode caCode caCode];
            end
            

            %--- Perform various initializations ------------------------------

            % define initial code frequency basis of NCO
            codeFreq      = settings.codeFreqBasis;
            % define residual code phase (in chips)
            remCodePhase  = 0.0;
            % define carrier frequency which is used over whole tracking period
            carrFreq      = channel(channelNr).acquiredFreq;
            carrFreqBasis = channel(channelNr).acquiredFreq;
            % define residual carrier phase
            remCarrPhase  = 0.0;

            %code tracking loop parameters
            oldCodeNco   = 0.0;
            oldCodeError = 0.0;

            %carrier/Costas loop parameters
            oldCarrNco   = 0.0;
            oldCarrError = 0.0;
            
            rawSignal_vec = [];

            %=== Process the number of specified code periods =================
            for loopCnt =  1:codePeriods

    %% GUI update -------------------------------------------------------------
                % The GUI is updated every 50ms. This way Matlab GUI is still
                % responsive enough. At the same time Matlab is not occupied
                % all the time with GUI task.
                if (rem(loopCnt, 50) == 0)
                    try
                        waitbar(loopCnt/codePeriods, ...
                                hwb, ...
                                ['Tracking: Ch ', int2str(channelNr), ...
                                ' of ', int2str(settings.numberOfChannels), ...
                                '; PRN#', int2str(channel(channelNr).PRN), ...
                                '; Completed ',int2str(loopCnt), ...
                                ' of ', int2str(codePeriods), ' msec']);                       
                    catch
                        % The progress bar was closed. It is used as a signal
                        % to stop, "cancel" processing. Exit.
                        disp('Progress bar closed, exiting...');
                        return
                    end
                end

    %% Read next block of data ------------------------------------------------            
                % Find the size of a "block" or code period in whole samples

                % Update the phasestep based on code freq (variable) and
                % sampling frequency (fixed)
                codePhaseStep = codeFreq / settings.samplingFreq;

                blksize = ceil((settings.codeLength-remCodePhase) / codePhaseStep);


               if (dataAdaptCoeff==1)
                    % Read in the appropriate number of samples to process this
                    % interation 
                    [rawSignal, samplesRead] = fread(fid, ...
                        dataAdaptCoeff*blksize, settings.dataType);

                    rawSignal = rawSignal';
               end

                if (dataAdaptCoeff==2)
                    % Read in the appropriate number of samples to process this
                    % interation 
                    [rawSignal, samplesRead] = fread(fid, ...
                        dataAdaptCoeff*blksize, settings.dataType);

                    rawSignal = rawSignal';
                    rawSignal1=rawSignal(1:2:end);
                    rawSignal2=rawSignal(2:2:end);
                    rawSignal = rawSignal1 + i .* rawSignal2;  %transpose vector
                end
                if (dataAdaptCoeff==4)   
                    [data, samplesRead] = fread(fid, [2, blksize], settings.dataType);
                    rawSignal = data(1,:) + data(2,:)*i; %Inphase and Quadrature
                    noise = sqrt(PowN/2).* (randn(1,length(data)) + 1j* randn(1,length(data)));
                    rawSignal = rawSignal + noise;
                    skip = 10;
                    rawSignal_vec = [rawSignal_vec, rawSignal];
                    if (loopCnt>skip) %
%                         index = randi((length(rawSignal_vec)-length(rawSignal)));
                        index = length(rawSignal_vec)-length(rawSignal)-randi(3);
                        rawSignal_total = rawSignal + 4.* rawSignal_vec(index:index+length(rawSignal)-1);
                    else
                        rawSignal_total = rawSignal;
                    end
                    rawSignal_vec = rawSignal;
                    rawSignal = rawSignal_total;
%                     cn0 = 37;
%                     rawSignal = rawSignal + sqrt(1/(10^(cn0/10)/settings.samplingFreq)).*(randn(1,blksize) + 1i.*randn(1,blksize));

                end
                 if (dataAdaptCoeff==8)    
                    [data, samplesRead] = fread(fid, [2, blksize], settings.dataType);
                    rawSignal = data(1,:) + data(2,:).*i; %Inphase and Quadrature

                 end

                % If did not read in enough samples, then could be out of 
                % data - better exit

                switch dataAdaptCoeff

                    case 4 
                    if (samplesRead ~= 2*blksize)
                        disp('Not able to read the specified number of samples  for tracking, exiting!')
                        fclose(fid);
                        return
                    end
                    case 8
                    if (samplesRead ~= 2*blksize)
                        disp('Not able to read the specified number of samples  for tracking, exiting!')
                        fclose(fid);
                        return
                    end
                    otherwise
                    if (samplesRead ~= dataAdaptCoeff*blksize)
                        disp('Not able to read the specified number of samples  for tracking, exiting!')
                        fclose(fid);
                        return
                    end
                end


    %% DLL discriminator
    Discriminator;
    %% Find PLL error and update carrier NCO ----------------------------------

                % Implement carrier loop discriminator (phase detector)
                carrError = atan(Q_P / I_P) / (2.0 * pi);

                % Implement carrier loop filter and generate NCO command
                carrNco = oldCarrNco + (tau2carr/tau1carr) * ...
                    (carrError - oldCarrError) + carrError * (PDIcarr/tau1carr);
                oldCarrNco   = carrNco;
                oldCarrError = carrError;

                % Modify carrier freq based on NCO command
                carrFreq = carrFreqBasis - carrNco;

                trackResults(channelNr).carrFreq(loopCnt) = carrFreq;

    %% Find DLL error and update code NCO -------------------------------------
                codeError = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / ...
                    (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));
                
                if(abs(codeError)>0.1)
                    1;
                end

                % Implement code loop filter and generate NCO command
                codeNco = oldCodeNco + (tau2code/tau1code) * ...
                    (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
                oldCodeNco   = codeNco;
                oldCodeError = codeError;

                % Modify code freq based on NCO command
                codeFreq = settings.codeFreqBasis - codeNco;

                trackResults(channelNr).codeFreq(loopCnt) = codeFreq;

    %% Record various measures to show in postprocessing ----------------------
                % Record sample number (based on 8bit samples)
               trackResults(channelNr).absoluteSample(loopCnt) =(ftell(fid))/dataAdaptCoeff;

               PRN_start_sample = PRN_start_sample+codeError; %javi
               trackResults(channelNr).PRN_start_sample(loopCnt)       = PRN_start_sample; %javi

                trackResults(channelNr).dllDiscr(loopCnt)       = codeError;
                trackResults(channelNr).dllDiscrFilt(loopCnt)   = codeNco;
                trackResults(channelNr).pllDiscr(loopCnt)       = carrError;
                trackResults(channelNr).pllDiscrFilt(loopCnt)   = carrNco;

                trackResults(channelNr).I_E(loopCnt) = I_E;
                trackResults(channelNr).I_P(loopCnt) = I_P;
                trackResults(channelNr).I_L(loopCnt) = I_L;
                trackResults(channelNr).Q_E(loopCnt) = Q_E;
                trackResults(channelNr).Q_P(loopCnt) = Q_P;
                trackResults(channelNr).Q_L(loopCnt) = Q_L;
            end % for loopCnt

           
            % If we got so far, this means that the tracking was successful
            % Now we only copy status, but it can be update by a lock detector
            % if implemented
            trackResults(channelNr).status  = channel(channelNr).status;        

        end % if a PRN is assigned
    end % for channelNr 
   

  
% Close the waitbar
close(hwb)
end