% clear
% close all
% clc

%% load parameter setting
addpath include             % The software receiver functions
addpath geoFunctions        % Position calculation related functions 
settings = initSettings();
% dataAdaptCoeff = 2;
numSV = 4;
PRN = [10 27 8 16];
% codePhase = zeros(1, 5);
% samplesPerCode = settings.samplingFreq * 1e-3;
% load("loadData/eph.mat")
% load("loadData/satInfo_navFile_PLAN_noSag.mat")
% load("caCode.mat")
satInfoRaw = load("loadData/satInfo_1ms.mat");
load("loadData/caCodesTable_25375.mat")
load("loadData/trackingResults_1ms.mat")
load("loadData/navSolutions_1ms.mat")
load("loadData/acqResults_1ms.mat")
load("loadData/pseudorage_correction.mat")
dtEst_1ms = load("loadData/dtEst_1ms.mat");

%% read raw data according to the solution time stamp
% leapSec = 18;
% startTransmitTime = leapSec + navSolutions.transmitTime(1) - ...
%     navSolutions.subFrameStart / 1000;
numSol = 75;%size(satInfo_navFile, 3);
readDataFlag = 0;
jumpMs = 797 - 10;
frame = navSolutions.subFrameStart;
frame = frame - min(frame) + 1;
if readDataFlag
    rawSignalDPE = cell(numSV,1);
    for iSv = 1:numSV
        skipNumberOfBytes = (jumpMs + frame(iSv) - 1) * 25375;
        [fid, message] = fopen(settings.fileName, 'rb');
        if (settings.fileType==2) 
            dataAdaptCoeff=2;
        end
        fseek(fid, dataAdaptCoeff * skipNumberOfBytes, 'bof');
    
        rawSignalDPE{iSv} = readDataDPE_20ms(fid, numSol, settings);
    end
end



truePosECEF = [4777973.177, 176346.307, 4207663.62];
% truePosECEF = [navSolutions.X(2), navSolutions.Y(2), navSolutions.Z(2)];
% trueCb = -8.4002e4 / settings.c;
trueCb = -4.28447e-4;%;-4.25743e-4;%;
trueVelECEF = [0 0 0];
trueCd = 0;
trueUsrInfo = [truePosECEF, trueCb, trueVelECEF, trueCd];

% subFrameStart = navSolutions.subFrameStart; 
Doppler = zeros(numSV, 80000);
for kSV = 1: numSV
    Doppler(kSV,:) = trackResults(kSV).carrFreq;
end

%% pre-set memory
posDPEECEF = zeros(numSol, 3);
dtEst = zeros(numSol,1);
%% load and extract raw data according to the tracking stage

% hwb0 = waitbar(0,'Solving position...');
% NSol = min(subFrameStart): 1000: min(subFrameStart) + 1000*(numSol-1);
for iSol = 1:numSol  
    % if mod(nSol,10) == 0
    %     waitbar(nSol/numSol, hwb0,  ...
    %         'Solving position...');
    % end
     %% extract raw signal
    jumpSol = 20;
    rawSignal = zeros(numSV, 25375*20);
    for iSv = 1:numSV
        rawSignal(iSv, :) = rawSignalDPE{iSv}(iSol,:);
    end
    % satInfo_nav = satInfo_navFile(:,:,iSol + jumpSol);
    % satInfo = satInfo_nav(:,1:7)';
    satInfo = satInfoRaw.satInfo(1:7,:,iSol);
    corr = correction(:,:,iSol + jumpSol);
    %% PDPE (grid search)
    % % generate search grid
    gridRangeCb = 1e-7;
    gridDensityCb = 1e-9;
    % centerPointCb = trueUsrInfo(4);
    if iSol == 1
        centerPointCb = trueUsrInfo(4);
    else
        centerPointCb = dtEst(iSol - 1) + 3.6e-8;
    end
    dtCan = (centerPointCb - gridRangeCb:gridDensityCb:centerPointCb + gridRangeCb);
    ndt = length(dtCan);
    rdt = zeros(1, ndt);

    %% generate local and perform correlation
    rMax_dt = zeros(ndt, 1);
    pos_dt = zeros(ndt, 3);
    for idt = 1: ndt
        usrInfo = [trueUsrInfo(1:3), dtCan(idt),...
                     trueUsrInfo(5:end), settings.IF];
        rdt(idt) = generateCAFs_20ms_corr(usrInfo, [satInfo; PRN], rawSignal,...
            caCodesTable, Doppler(:,(iSol-1)*1000 + 1), corr);
    end
    [~, maxInd] = max(rdt);
    % posDPEECEF(iSol,:) = pos_dt(maxInd, :);
    dtEst(iSol) = dtCan(maxInd);
    %% PDPE (generalized pattern search)
    % usrInfo = [trueUsrInfo(1) - 100, trueUsrInfo(2) + 100, trueUsrInfo(3:end), settings.IF];
    % % usrInfo = [trueUsrInfo, settings.IF];
    % satInfoPRN = [satInfo(:,:,nSol); PRN];
    % searchD = [eye(2), [-1; -1]];
    % stepSize = 20;
    % stepDecay = 0.6;
    % convergence = 0.01;
    % [localMax, cost, searchPath] = generalizedPatternSearch2d(usrInfo, satInfoPRN, rawSignal, caCodesTable,...
    % timeDiff, searchD, stepSize, stepDecay, convergence);
    % posDPEECEF(nSol, :) = [localMax(1:2), truePosECEF(3)];

    %% visualize CAF (2D)
    figure(102), plot(dtCan, rdt)
    % xlim([-4.285e-4, -4.257e-4])
    text(dtCan(maxInd),rdt(maxInd),sprintf('%d',iSol))
    hold on
    
end
% close(hwb0)
%% plot solutions
save("result/dtEst.mat", "dtEst")




