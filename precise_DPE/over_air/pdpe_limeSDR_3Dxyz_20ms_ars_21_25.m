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
dtEst_20ms = load("loadData/dtEst_20ms.mat");

%% read raw data according to the solution time stamp
% leapSec = 18;
% startTransmitTime = leapSec + navSolutions.transmitTime(1) - ...
%     navSolutions.subFrameStart / 1000;
numSol = 75;%size(satInfo_navFile, 3);
% readDataFlag = 0;
% jumpMs = 797 - 10;
% frame = navSolutions.subFrameStart;
% frame = frame - min(frame) + 1;
% if readDataFlag
%     rawSignalDPE = cell(numSV,1);
%     for iSv = 1:numSV
%         skipNumberOfBytes = (jumpMs + frame(iSv) - 1) * 25375;
%         [fid, message] = fopen(settings.fileName, 'rb');
%         if (settings.fileType==2) 
%             dataAdaptCoeff=2;
%         end
%         fseek(fid, dataAdaptCoeff * skipNumberOfBytes, 'bof');
% 
%         rawSignalDPE{iSv} = readDataDPE_20ms(fid, numSol, settings);
%     end
% end

% record raw signal
% for iSol = 1:numSol
%     rawSignal = zeros(numSV, 25375*20);
%     for iSv = 1:numSV
%         rawSignal(iSv, :) = rawSignalDPE{iSv}(iSol,:);
%     end
%     fileName = sprintf('loadData/rawSignal/rawSignal_sol%d.mat',iSol);
%     save(fileName, "rawSignal")
% end


%% true information
truePosECEF = [4777973.177, 176346.307, 4207663.62];
% truePosECEF = [[4777962.17700000	176358.307000000	4207657.62000000]];
% truePosECEF = [4777982.17700000	176346.307000000	4207662.62000000];
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

%% setting for generalized pattern search
dmax = 1;
dmin = 0.01;
contraction = 2;
nIter = 300;

curFile = mfilename;
sol1 = round(str2double(curFile(end - 4:end - 3)));
sol2 = round(str2double(curFile(end - 1:end)));

for iSol = sol1:sol2
    iSol
    % if mod(nSol,10) == 0
    %     waitbar(nSol/numSol, hwb0,  ...
    %         'Solving position...');
    % end
    %% extract raw signal
    jumpSol = 20;
    % rawSignal = zeros(numSV, 25375*20);
    % for iSv = 1:numSV
    %     rawSignal(iSv, :) = rawSignalDPE{iSv}(iSol,:);
    % end
    fileName = sprintf('loadData/rawSignal/rawSignal_sol%d.mat',iSol);
    signal = load(fileName);
    rawSignal = signal.rawSignal;
    satInfo = satInfoRaw.satInfo(1:7,:,iSol);
    corr = correction(:,:,iSol + jumpSol);
    %% DPE (ARS)
    dt = dtEst_20ms.dtEst(iSol);
    % generate z candidate
    zRange = 0.12;
    zDensity = 0.01;
    zCan = (truePosECEF(3) - zRange + zDensity: zDensity: truePosECEF(3) + zRange);
    % zCan = grid1d(truePosECEF(3), zRange, zDensity);
    nz = length(zCan);
    % hwb = waitbar(0, 'Geps for candidates z');
    %% generate local and perform correlation
    rMax_z = zeros(nz, 1);
    pos_z = zeros(nz, 2);
    parfor iz = 1: nz
        iz
        % waitbar(iz/nz, hwb, sprintf('Geps for candidates z %d/%d',iz, nz));
        usrInfo = [trueUsrInfo(1:2), zCan(iz), dt,...
                     trueUsrInfo(5:end), settings.IF];
        [pos_z(iz,:), rMax_z(iz)] = ars2d_pdpe_20ms(usrInfo, [satInfo; PRN], rawSignal, caCodesTable,...
                        Doppler, corr, dmax, dmin, contraction, nIter);

    end
    [~, maxInd] = max(rMax_z);
    posDPEECEF(iSol,:) = [pos_z(maxInd, :), zCan(maxInd)];
    % close(hwb)

end
% save("dtEst_1ms.mat", "dtEst")
% close(hwb0)
%% plot solutions
% figure, plot(dtCan, rdt)
% [maxR, maxInd] = max(r);
% maxDt = dtCan(maxInd);
% save("result/rMax2d.mat", "rMax2d")
% save("result/dtEst_1_2.mat", "dtEst")
saveName = sprintf("result/posDPEECEF_20ms_ars_%d_%d.mat",sol1, sol2);
save("saveName", "posDPEECEF")
% solPlot(settings, truePosECEF, navSolutions, posDPEECEF, 53)
% settings = initSettings();
% plotNavigation(navSolutions, settings,1)



