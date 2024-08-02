clear
clc
close all
%%
addpath (genpath ('goGPS'));

roverObservationFile = 'SABA1400.20O';
useGps = 1;   % Use frequencies L1, L2, L5
useGal = 1;   % Use frequencies E1, E5a, E5b
[constellations] = goGNSS.initConstellation(useGps,0,useGal,0,0,0);
[codeL1Rover, phaseL1Rover, codeL2Rover, phaseL2Rover, dopplerL1Rover,...
    dopplerL2Rover, snrL1Rover, snrL2Rover, timeGpsRefRover, timeRover, weekRover, dateRover,...
    posReceiverRover, obsIntervalRover, ~, ~, C1RoverIsUsed] = load_RINEX_obs(roverObservationFile, constellations);
