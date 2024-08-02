% Script postProcessing.m processes the raw signal from the specified data
% file (in settings) operating on blocks of 37 seconds of data.
%
% First it runs acquisition code identifying the satellites in the file,
% then the code and carrier for each of the satellites are tracked, storing
% the 1msec accumulations.  After processing all satellites in the 37 sec
% data block, then postNavigation is called. It calculates pseudoranges
% and attempts a position solutions. At the end plots are made for that
% block of data.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis, Dennis M. Akos
% Some ideas by Dennis M. Akos
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

%                         THE SCRIPT "RECIPE"
%
% The purpose of this script is to combine all parts of the software
% receiver.
%
% 1.1) Open the data file for the processing and seek to desired point.
%
% 2.1) Acquire satellites
%
% 3.1) Initialize channels (preRun.m).
% 3.2) Pass the channel structure and the file identifier to the tracking
% function. It will read and process the data. The tracking results are
% stored in the trackResults structure. The results can be accessed this
% way (the results are stored each millisecond):
% trackResults(channelNumber).XXX(fromMillisecond : toMillisecond), where
% XXX is a field name of the result (e.g. I_P, codePhase etc.)
%
% 4) Pass tracking results to the navigation solution function. It will
% decode navigation messages, find satellite positions, measure
% pseudoranges and find receiver position.
%
% 5) Plot the results.

%% Initialization =========================================================
disp ('Starting processing...');

[fid, message] = fopen(settings.fileName, 'rb');

%Initialize the multiplier to adjust for the data type
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

%If success, then process the data
if (fid > 0)
    
    % Move the starting point of processing. Can be used to start the
    % signal processing at any point in the data record (e.g. good for long
    % records or for signal processing in blocks).
    fseek(fid, dataAdaptCoeff*settings.skipNumberOfBytes, 'bof'); 
    % fseek(fid, 493699694, 'bof'); 
%% Acquisition ============================================================

    % Do acquisition if it is not disabled in settings or if the variable
    % acqResults does not exist.
    if ((settings.skipAcquisition == 0) || ~exist('acqResults', 'var'))
        
        % Find number of samples per spreading code
        samplesPerCode = round(settings.samplingFreq / ...
                           (settings.codeFreqBasis / settings.codeLength));
        
        % Read data for acquisition. 11ms of signal are needed for the fine
        % frequency estimation
        if (dataAdaptCoeff==1)         
             data  = fread(fid, dataAdaptCoeff*22*samplesPerCode, settings.dataType)';
        end 
        if (dataAdaptCoeff==2)    
            [data, count] = fread(fid, [2, 22*samplesPerCode], settings.dataType);
            data = data(1,:) + data(2,:).*1i; %Inphase and Quadrature
        end

        if (dataAdaptCoeff==4)    
            [data, count] = fread(fid, [2, 22*samplesPerCode], settings.dataType);
            data = data(1,:) + data(2,:)*1i; %Inphase and Quadrature
        end

        if (dataAdaptCoeff==8)    
            [data, count] = fread(fid, [2, 22*samplesPerCode], settings.dataType);
            data = data(1,:) + data(2,:).*1i; %Inphase and Quadrature
        end
        
        %--- Do the acquisition -------------------------------------------
        disp ('   Acquiring satellites...');
        if (settings.highsensitivity_acquisition==1)
           acqResults = acquisition_10ms(data, settings); %Javier: High Sensitivity acquisition
        else
            acqResults = acquisition(data, settings);
        end
        plotAcquisition(acqResults);
    end

%% Initialize channels and prepare for the run ============================

    % Start further processing only if a GNSS signal was acquired (the
    % field FREQUENCY will be set to 0 for all not acquired signals)
    if (any(acqResults.carrFreq))
        channel = preRun(acqResults, settings);
        showChannelStatus(channel, settings);
    else
        % No satellites to track, exit
        disp('No GNSS signals detected, signal processing finished.');
        trackResults = [];
        return;
    end

%% Track the signal =======================================================
    startTime = now;
    disp (['   Tracking started at ', datestr(startTime)]);

    % Process all channels for given data block
    % do a normal tracking to get the subframe start
    % [trackResults, channel] = tracking(fid, channel, settings); 
    % do another tracking to record the raw data for DPE
    [trackResults, channel] = tracking_20ms(fid, channel, settings);

    % % Close the data file
    % fclose(fid);
    
    disp(['   Tracking is over (elapsed time ', ...
                                        datestr(now - startTime, 13), ')'])     

    % Auto save the acquisition & tracking results to a file to allow
    % running the positioning solution afterwards.
    disp('   Saving Acq & Tracking results to file "trackingResults.mat"')
    save('trackingResults', ...
                      'trackResults', 'settings', 'acqResults', 'channel'); 

    save("trackingResults_1ms.mat", "trackResults")
%% Calculate navigation solutions =========================================
    disp('   Calculating navigation solutions...');
    [navSolutions, eph] = postNavigation_20ms(trackResults, settings, acqResults);
    % [navSolutions, ~, subFrameStart] = postNavigation_record(trackResults, settings);
    % record for PDPE
    save("eph.mat", "eph")
    save("navSolutions_1ms.mat", "navSolutions")
    disp('   Processing is complete for this data block');
    
%% Record the raw tracking satellite signal based on the preamble
%     [rawSignal_record, recordFlag] = tracking_record(fid, channel, settings, subFrameStart);
%     % Close the data file
%     fclose(fid);
%     % Record the time stamps, satellite positions and raw signals
%     RawSignalAndNavi.recordFlag = recordFlag.';
% %     RawSignalAndNavi.chipPreamble = chipPreamble;
%     RawSignalAndNavi.numSol = length(recordFlag);
%     numSat = size(navSolutions.channel.PRN, 1);
%     RawSignalAndNavi.numSat = numSat;
% %     diffFrame = subFrameStart - min(subFrameStart);
% %     RawSignalAndNavi.delayedTime = repmat(navSolutions.transmitTime.', 1, numSat) + ...
% %         repmat(diffFrame, length(navSolutions.transmitTime), 1);
%     RawSignalAndNavi.satPos = navSolutions.satPositions;
%     RawSignalAndNavi.satVel = navSolutions.satVelocities;
%     RawSignalAndNavi.satClkCorr = navSolutions.satClkCorr;
%     RawSignalAndNavi.rawSignal = rawSignal_record;
%     RawSignalAndNavi.PRN = navSolutions.channel.PRN.';
%     RawSignalAndNavi.posEst = [navSolutions.X.', navSolutions.Y.', navSolutions.Z.'];
%     RawSignalAndNavi.cbEst = navSolutions.dt.';
% 
%     RawSignalAndNavi.delayedTime = navSolutions.delayedTime;
%     RawSignalAndNavi.satPositionsDelayed = navSolutions.satPositionsDelayed;
%     RawSignalAndNavi.satVelocitiesDelayed = navSolutions.satVelocitiesDelayed;
%     RawSignalAndNavi.satClkCorrDelayed = navSolutions.satClkCorrDelayed;
%     for nChn = 1: numSat
%         RawSignalAndNavi.Dp(:, nChn) = trackResults(nChn).carrFreq(subFrameStart(nChn): ...
%             10:(subFrameStart(nChn) + 10*(RawSignalAndNavi.numSol-1))); 
%     end
%     RawSignalAndNavi.travelTime = navSolutions.channel.correctedP'/ settings.c * 1e3;
% 
%     disp('   Saving RawSignal and Navi Solutions to file "RawSignalAndNavi.mat"')
%     save('RawSignalAndNavi_SPAN.mat',"RawSignalAndNavi");
%% Plot all results ===================================================
    % disp ('   Ploting results...');
    % if settings.plotTracking
    %     plotTracking(1:settings.numberOfChannels, trackResults, settings);
    % end

    % plot RAW pseudoranges
%     figure;
%     plot(navSolutions.transmitTime,navSolutions.channel(:).rawP.')
%     title('RAW Pseudoranges vs. TOW');
%     ylabel('Pseudorange [m]');
%     xlabel('GPS TOW [s]');

    plotNavigation(navSolutions, settings,1);
    writematrix(navSolutions.solECEF,'loadData/navSolECEF.txt','Delimiter','tab');
    % GoogleEarth(navSolutions,settings)

    disp('Post processing of the signal is over.');

else
    % Error while opening the data file.
    error('Unable to read file %s: %s.', settings.fileName, message);
end % if (fid > 0)
