function [corr, codeL1Rover] = generate_correction(obsFile, satInfo, PRN, settings)
%GENERATE_CORRECTION generate corrections from precise obs file
% corr: nSat by 4 by nEpoch
%       col 1: pseudorange correction in meter (including receiver's clock bias)
%       col 2: carrier phase correction in cycle
%       col 3: transmitted time in GPS TOW
%       col 4: sat PRN
addpath (genpath ('include/goGPS'));

truePos = settings.truePosition.ECEF;
trueVel = 0;
useGps = 1;   % Use frequencies L1, L2, L5
useGal = 1;   % Use frequencies E1, E5a, E5b
[constellations] = goGNSS.initConstellation(useGps,0,useGal,0,0,0);
[codeL1Rover, phaseL1Rover, ~, ~, dopplerL1Rover,...
    ~, snrL1Rover, ~, timeGpsRefRover, timeRover, weekRover, dateRover,...
    posReceiverRover, ~, ~, ~, ~] = load_RINEX_obs(obsFile, constellations);
TOW = mod(timeRover + 18, weekRover*7*24*60*60);

nEpoch = size(satInfo,3);
nSat = length(PRN);
corr = zeros(nSat, 4, nEpoch);
lambda = settings.c/settings.GPS_L1;
for iEpoch = 1: nEpoch
    satPositions = satInfo(:,1:3,iEpoch);
    satVelocities = satInfo(:,4:6,iEpoch);
    satClkCorr = satInfo(:,7,iEpoch);

    pseudoranges = vecnorm(satPositions - repmat(truePos, nSat, 1), 2, 2) - settings.c * satClkCorr;
    % unitvector = (satPositions' - truePos) ./ pseudoranges;
    carrierphase = pseudoranges/ lambda;

    pseudorangesObs = codeL1Rover(PRN, iEpoch);
    carrierphaseObs = phaseL1Rover(PRN, iEpoch);

    pseudorangesCorr = pseudoranges - pseudorangesObs;
    carrierphaseCorr = carrierphase - carrierphaseObs;
    corr(:,1,iEpoch) = pseudorangesCorr;
    corr(:,2,iEpoch) = carrierphaseCorr;
    corr(:,3,iEpoch) = TOW(iEpoch) * ones(4,1);
    corr(:,4,iEpoch) = PRN;

end




