clear
close all
clc

%% load parameter setting
addpath include             % The software receiver functions
addpath geoFunctions        % Position calculation related functions 
settings = initSettings();
% dataAdaptCoeff = 2;
numSV = 4;
PRN = [10 27 8 16];
coherent = 10;
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
jumpMs = 797;
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
    
        rawSignalDPE{iSv} = readDataDPE_coherent(fid, numSol, settings, coherent);
    end
end


%% true information
truePosECEF = [4777973.177, 176346.307, 4207663.62];
% truePosECEF = [navSolutions.X(2), navSolutions.Y(2), navSolutions.Z(2)];
% trueCb = -8.4002e4 / settings.c;
trueCb = -4.28447e-4;%-4.25743e-4;%;
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
for iSol = 1:2
    % if mod(nSol,10) == 0
    %     waitbar(nSol/numSol, hwb0,  ...
    %         'Solving position...');
    % end
    %% extract raw signal
    jumpSol = 20;
    rawSignal = zeros(numSV, 25375*coherent);
    for iSv = 1:numSV
        rawSignal(iSv, :) = rawSignalDPE{iSv}(iSol,:);
    end
    % satInfo_nav = satInfo_navFile(:,:,iSol + jumpSol);
    % satInfo = satInfo_nav(:,1:7)';
    satInfo = satInfoRaw.satInfo(1:7,:,iSol);
    corr = correction(:,:,iSol + jumpSol);
    %% PDPE (grid search)
    % % generate search gridz
    gridRangeCb = 1e-8;
    gridDensityCb = 1e-9;
    % centerPointCb = trueUsrInfo(4);
    % if iSol == 1
    %     centerPointCb = trueUsrInfo(4);
    % else
    %     centerPointCb = dtEst(iSol - 1) + 3.6e-8;
    % end
    centerPointCb = dtEst_1ms.dtEst(iSol);
    dtCan = (centerPointCb - gridRangeCb:gridDensityCb:centerPointCb + gridRangeCb);
    

    gridRangeXy = 30;
    gridDensityXy = 1;
    centerPointXyz = trueUsrInfo;
    [xCan, yCan, zCan] = grid3d(centerPointXyz, gridRangeXy, gridDensityXy);

    nx = length(xCan);
    ny = length(yCan);
    nz = length(zCan);
    ndt = length(dtCan);
    r = zeros(nx,ny);
    rdt = zeros(1, ndt);
    rMax = 0;
    rMax2d = zeros(nx, ny, numSol);
    % hwb = waitbar(0,'Correlation...');
    %% generate local and perform correlation
    rMax_dt = zeros(ndt, 1);
    pos_dt = zeros(ndt, 3);
    for idt = 1: ndt
        % if mod(idt, 100) == 0
        %     waitbar(idt/ndt, hwb, 'Correlation...');
        % end
        rMax_zsave = 0;
        for iz = 1:nz
            for ix = 1: nx
                parfor iy = 1: ny                
                    usrInfo = [xCan(ix), yCan(iy), zCan(iz), dtCan(idt),...
                        trueUsrInfo(5:end), settings.IF];
                    % usrInfo = [trueUsrInfo, settings.IF];
                    r(ix, iy) = generateCAFs_coherent_corr(usrInfo, [satInfo; PRN], rawSignal,...
                        caCodesTable, Doppler(:,(iSol-1)*1000 + 1), corr, coherent);
                end
            end
            rMax_z = max(r(:));
            if rMax_z > rMax_zsave
                rMax_zsave = rMax_z;
                [mMax,nMax] = find(r == rMax_z);
                pos_dt(idt,:) = [xCan(nMax), yCan(mMax), zCan(iz)];
                rMax_dt(idt) = rMax_z;
            end    
        end
        % usrInfo = [trueUsrInfo(1:3), dtCan(idt),...
        %              trueUsrInfo(5:end), settings.IF];
        % rdt(idt) = generateCAFs_1ms_corr(usrInfo, [satInfo; PRN], rawSignal,...
        %     caCodesTable, Doppler(:,(iSol-1)*1000 + 1), corr);
    end
    [~, maxInd] = max(rMax_dt);
    posDPEECEF(iSol,:) = pos_dt(maxInd, :);
    dtEst(iSol) = dtCan(maxInd);

    % close(hwb)
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
    % figure(102), plot(dtCan, rdt)
    % xlim([-4.285e-4, -4.257e-4])
    % hold on
    % r = rMax2d;
    % [mMax,nMax] = find(r == max(r(:)));
    % plot2DCAF(xCan, yCan, r, truePosECEF);
    % posDPEECEF(nSol, :) = [xCan(nMax), yCan(mMax), truePosECEF(3)];
    % figure(100), hold on;
    % plot3(searchPath(:,1), searchPath(:,2), searchPath(:,3)+ 1e9,...
    %     "Color", "#A2142F", "LineWidth",1.5)
    % [maxR_dt, maxInd_dt] = max(rdt);
    % dtEst(iSol) = dtCan(maxInd_dt);

    
end
% save("dtEst_1ms.mat", "dtEst")
% close(hwb0)
%% plot solutions
% figure, plot(dtCan, rdt)
% [maxR, maxInd] = max(r);
% maxDt = dtCan(maxInd);
% save("result/rMax2d.mat", "rMax2d")
save("result/dtEst_1_2.mat", "dtEst")
save("result/posDPEECEF_1_2.mat", "posDPEECEF")
pos2SPECEF = [navSolutions.X',navSolutions.Y',navSolutions.Z'];
% solPlot(truePosECEF, pos2SPECEF(1:5,:), posDPEECEF(1:5,:))





