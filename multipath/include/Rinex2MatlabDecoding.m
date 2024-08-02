function [ SatelliteInformation, obsIntervalRover, nEpochs ] = Rinex2MatlabDecoding( varargin )
% -------------------------------------------------------------------------
%  DLR Institute of Communications and Navigation
%
%  Authors: Daniel Medina, Carlos Alvarez Merino
%
%                       Rinex2MatlabDecoding.m
%
%  BLABLABLAthisfunctiondoesthis....
%
% DEPENDENCES
%   - goGPS software for reading RINEX files, computing satellite
%   positions and estimate ionospheric and tropospheric corrections.
%
% REFERENCES
%   - goGPS webpage: https://github.com/goGPS-Project/goGPS_MATLAB
%
% TO DO LIST:
% - Add the simulation for the RTK (it is just a copy-paste from a already
% existing code, but it needs some cleaning up)
% 
% Date: 25 July 2019
% -------------------------------------------------------------------------

%% Handle input arguments
% Default values
useSimulation = 1;
useGps = 1;
useGal = 0;
elevationMask = 15;
gpsUnhealthy = [];
galUnhealthy = [14 18];
if nargin > 2
    for idx = 1:2:length(varargin)
        if ischar( varargin{idx} )
            switch varargin{idx}
                case 'roverObservationFile'
                    roverObservationFile = varargin{idx+1};
                case 'baseObservationFile'
                    baseObservationFile = varargin{idx+1}; % Yes, base and rover are exchanged. This is in purpose at the moment, but it will change soon
                case 'useGps'
                    useGps = varargin{idx+1};
                case 'useGal'
                    useGal = varargin{idx+1};
                case 'navigationFile'
                    navigationFile = varargin{idx+1};
                case 'nEpochs'
                    nEpochs2Use= varargin{idx+1};
                case 'useSimulation'
                    useSimulation= varargin{idx+1};
                case 'elevationMask'
                    elevationMask= varargin{idx+1};
                case 'gpsUnhealthy'
                    gpsUnhealthy= varargin{idx+1};
                case 'galUnhealthy'
                    galUnhealthy= varargin{idx+1};
                case 'nEpochs2Use'
                    nEpochs2Use= varargin{idx+1};
                otherwise
                    error(['Unrecognized variable: ' varargin{idx}])
            end
        end
    end
elseif nargin<4
    error('Wrong number of inputs')
end

% Some error checking
if exist('baseObservationFile') == false
    error('Observation file is missing!')
end
if exist('navigationFile') == false 
    error('Navigation data file is missing!')
end
if useSimulation == true && exist('baseObservationFile') == false 
    error('The observation file of the base station is missing!')
end

% Arranging the GNSS frequencies to be used
if length( useGps ) == 1 
    useGpsL1 = useGps(1);
    useGpsL2 = 0;
    useGpsL5 = 0;
elseif length( useGps ) == 3
    useGpsL1 = useGps(1);
    useGpsL2 = useGps(2);
    useGpsL5 = useGps(3);
end
if length( useGal ) == 1 
    useGalE1 = useGal(1);
    useGalE5a = 0;
    useGalE5b = 0;
elseif length( useGal ) == 3
    useGalE1 = useGal(1);
    useGalE5a = useGal(2);
    useGalE5b = useGal(3);
end
frequencies2Use = [useGpsL1, useGpsL2, useGpsL5, useGalE1, useGalE5a, useGalE5b];
simulatedL5  = useGpsL5 && useSimulation;
simulatedE5a = useGalE5a && useSimulation;
simulatedE5b = useGalE5b && useSimulation;


%% GNSS Constants
SPEED_OF_LIGHT = 299792458;

WAVELENGTH_L1 = SPEED_OF_LIGHT/(goGNSS.FL1*1e6);
WAVELENGTH_L2 = SPEED_OF_LIGHT/(goGNSS.FL2*1e6);
WAVELENGTH_L5 = SPEED_OF_LIGHT/(goGNSS.FL5*1e6);
WAVELENGTH_E1 = SPEED_OF_LIGHT/(goGNSS.FE1*1e6);
WAVELENGTH_E5a = SPEED_OF_LIGHT/(goGNSS.FE5a*1e6);
WAVELENGTH_E5b = SPEED_OF_LIGHT/(goGNSS.FE5b*1e6);
WAVELENGTHS(1) = WAVELENGTH_L1; 
WAVELENGTHS(2) = WAVELENGTH_L2; 
WAVELENGTHS(3) = WAVELENGTH_L5;
WAVELENGTHS(4) = WAVELENGTH_E1; 
WAVELENGTHS(5) = WAVELENGTH_E5a; 
WAVELENGTHS(6) = WAVELENGTH_E5b;

nFreq = length(WAVELENGTHS);


%% 1) Load GNSS data
useGps = max(useGps);
useGal = max(useGal);
[constellations] = goGNSS.initConstellation(useGps,0,useGal,0,0,0);
nSatTot = constellations.nEnabledSat;
n_sys = useGps + useGal; % Number of constellations used

% Load the observation file for the Rover station
[codeL1Rover, phaseL1Rover, codeL2Rover, phaseL2Rover, dopplerL1Rover, dopplerL2Rover, snrL1Rover, snrL2Rover, timeGpsRefRover, timeRover, weekRover, dateRover, posReceiverRover, obsIntervalRover, ~, ~, C1RoverIsUsed] = load_RINEX_obs(roverObservationFile, constellations);

[Eph, iono] = load_RINEX_nav(navigationFile, constellations, 0);
%Rover Station
% obsParametersRover is a combination of C1, L1, D1, S1, C2, L2, D2, S2, C5, L5, D5, S5
obsParametersRover.C1 = codeL1Rover;
obsParametersRover.L1 = phaseL1Rover;
obsParametersRover.D1 = dopplerL1Rover;
obsParametersRover.S1 = snrL1Rover;
obsParametersRover.C2 = codeL2Rover;
obsParametersRover.L2 = phaseL2Rover;
obsParametersRover.D2 = dopplerL2Rover;
obsParametersRover.S2 = snrL2Rover;
%We suppose that L5, C5 and D5 will be implemented in the future
obsParametersRover.C5 = zeros(size(codeL1Rover));
obsParametersRover.L5 = zeros(size(phaseL1Rover));
obsParametersRover.D5 = zeros(size(dopplerL1Rover));
obsParametersRover.S5 = zeros(size(snrL1Rover));

nEpochsRover = length(timeRover);
[phiRover, lamRover, hRover] = cart2geod(posReceiverRover(1), posReceiverRover(2), posReceiverRover(3));
phiRover = phiRover * 180 / pi;
lamRover = lamRover * 180 / pi;
dtB = zeros(nEpochsRover,1); % receiver clock error

% Load the observation file for the base (only if using real data)
if useSimulation == false
    
    % Load the observation file for the rover
    [codeL1Base, phaseL1Base, codeL2Base, phaseL2Base, dopplerL1Base, dopplerL2Base, snrL1Base, snrL2Base, timeGpsRefBase, timeBase, weekBase, dateBase, posReceiverBase, obsIntervalBase, ~, ~, C1BaseIsUsed] = load_RINEX_obs(baseObservationFile, constellations);
    
    %Base
    % obsParametersBase is a combination of C1, L1, D1, S1, C2, L2, D2, S2, C5, L5, D5, S5
    obsParametersBase.C1 = codeL1Base;
    obsParametersBase.L1 = phaseL1Base;
    obsParametersBase.D1 = dopplerL1Base;
    obsParametersBase.S1 = snrL1Base;
    obsParametersBase.C2 = codeL2Base;
    obsParametersBase.L2 = phaseL2Base;
    obsParametersBase.D2 = dopplerL2Base;
    obsParametersBase.S2 = snrL2Base;
    %We suppose that L5, C5 and D5 will be implemented in the future
    obsParametersBase.C5 = zeros(size(codeL1Base));
    obsParametersBase.L5 = zeros(size(phaseL1Base));
    obsParametersBase.D5 = zeros(size(dopplerL1Base));
    obsParametersBase.S5 = zeros(size(snrL1Base));
    
    nEpochsBase = length(timeBase);
    [phiBase, lamBase, hBase] = cart2geod(posReceiverBase(1), posReceiverBase(2), posReceiverBase(3));
    phiBase = phiBase * 180 / pi;
    lamBase = lamBase * 180 / pi;
    dtR = zeros(nEpochsBase,1); % receiver clock error
end


if exist('nEpochs2Use','var') == false
    nEpochs2Use = nEpochsRover;
end
if exist('nEpochsBase','var') == true
    nEpochs = min( [ nEpochsRover, nEpochs2Use, nEpochsBase ] );
else
    nEpochs = min( [ nEpochsRover, nEpochs2Use ] );
end



%% Loading Satellite Ephemeris
nSatRover = zeros(nEpochs,1);
nSatBase = zeros(nEpochs,1);
typeFreqRover = NaN(length(WAVELENGTHS),nEpochs);
typeFreqBase = NaN(length(WAVELENGTHS),nEpochs);
iEpoch = 0;

while iEpoch < nEpochs
    
    iEpoch = iEpoch + 1;
          
    %%%%%%%%%%%%%%%%%% Base %%%%%%%%%%%%%%%%%%
    %Cell structure = [Gps-L1, Gps-L2, Gps-L5, GAL-E1, GAL-E5a, GAL-E5b]
    [ obsParametersRoverAux, posSatRover, velSatRover, ionoCorrectionRover, tropoCorrectionRover, clkSatRover, ...
        elevSatRover, distSatRover, nSatRover(iEpoch), labelSatRover, satPRN_R,  nSatConstRover, ...
        typeFreqRover(:,iEpoch), nSatConstExtendedRover, AzimuthRoverase, elevationMask] = ...
        SVPositionRINEXdecoder( iEpoch, Eph, obsParametersRover, timeGpsRefRover, nSatTot, posReceiverRover, ...
        phiRover, lamRover, hRover, iono, constellations, elevationMask );
    
    if simulatedL5 == true
        posSatRover{3} = posSatRover{1};    
        elevSatRover{3} = elevSatRover{1};
        typeFreqRover(3,iEpoch) = typeFreqRover(1,iEpoch);
    end
    if simulatedE5a == true
        posSatRover{5} = posSatRover{4};
        elevSatRover{5} = elevSatRover{4};
        typeFreqRover(5,iEpoch) = typeFreqRover(4,iEpoch);
    end
    if simulatedE5b == true  
        posSatRover{6} = posSatRover{4};
        elevSatRover{6} = elevSatRover{4};
        typeFreqRover(6,iEpoch) = typeFreqRover(4,iEpoch);
    end    
    
    for iFreq = 1:nFreq
        if frequencies2Use(iFreq) == false
            posSatRover{iFreq} = [];
            labelSatRover{iFreq} = [];
            elevSatRover{iFreq} = [];
            if useSimulation == false
                posSatRover{iFreq} = [];
                labelSatRover{iFreq} = [];
                elevSatRover{iFreq} = [];
            end
            typeFreqRover(iFreq,iEpoch) = 0;       
        end
    end
     
    for iFreq = 1:nFreq % Eliminate unhealthy satellites
        % GPS Satellites broken XX
        if iFreq < 4
            satellitesUnhealthy = gpsUnhealthy;
        else
            satellitesUnhealthy = galUnhealthy;
        end
           [~,indexSatBroken,~] = intersect(labelSatRover{iFreq},satellitesUnhealthy );
           if isempty(indexSatBroken) == false
               elevSatRover{iFreq}(indexSatBroken) = [];
               posSatRover{iFreq}(indexSatBroken,:) = [];
               labelSatRover{iFreq}(indexSatBroken) = [];
               typeFreqRover(iFreq,iEpoch) = typeFreqRover(iFreq,iEpoch) - length(indexSatBroken);
               obsParametersRoverAux.C{iFreq}(indexSatBroken) = [];
               obsParametersRoverAux.L{iFreq}(indexSatBroken) = [];
           end
    end % Eliminate unhealthy satellites

    SatelliteInformation(iEpoch).obsParametersRover = obsParametersRoverAux;
    SatelliteInformation(iEpoch).posSatRover = posSatRover;      %SV position
    SatelliteInformation(iEpoch).labelSatRover = labelSatRover;       %Number of the SV
    SatelliteInformation(iEpoch).elevSatRover = elevSatRover;     %Elev of the SV
      
    if useSimulation == false
        %%%%%%%%%%%%%%%%%% Rover %%%%%%%%%%%%%%%%%%
        %Cell structure = [Gps-L1, Gps-L2, Gps-L5, GAL-E1, GAL-E5a, GAL-E5b]
        [ obsParametersBaseAux, posSatBase, velSatBase, ionoCorrectionBase, tropoCorrectionBase, clkSatBase, ...
            elevSatBase, distSatBase, nSatBase(iEpoch), labelSatBase, satPRN_B,  nSatConstBase, ...
            typeFreqBase(:,iEpoch), nSatConstExtendedBase, AzimuthBase, elevationMask ] = ...
            SVPositionRINEXdecoder( iEpoch, Eph, obsParametersBase, timeGpsRefBase, nSatTot, posReceiverBase, ...
            phiBase, lamBase, hBase, iono, constellations, elevationMask );
        
        for iFreq = 1:nFreq % Eliminate unhealthy satellites
            % GPS Satellites broken XX
            if iFreq < 4
                satellitesUnhealthy = gpsUnhealthy;
            else
                satellitesUnhealthy = galUnhealthy;
            end
               [~,indexSatBroken,~] = intersect(labelSatBase{iFreq},satellitesUnhealthy );
               if isempty(indexSatBroken) == false
                   elevSatBase{iFreq}(indexSatBroken) = [];
                   posSatBase{iFreq}(indexSatBroken,:) = [];
                   labelSatBase{iFreq}(indexSatBroken) = [];
                   typeFreqBase(iFreq,iEpoch) = typeFreqBase(iFreq,iEpoch) - length(indexSatBroken);
                   obsParametersBaseAux.C{iFreq}(indexSatBroken) = [];
                   obsParametersBaseAux.L{iFreq}(indexSatBroken) = [];
               end
        end % Eliminate unhealthy satellites
        
        SatelliteInformation(iEpoch).obsParametersBase = obsParametersBaseAux;
        SatelliteInformation(iEpoch).posSatBase = posSatBase;      %SV position
        SatelliteInformation(iEpoch).labelSatBase = labelSatBase;       %Number of the SV 
        
        %Adjutsing just for the intersection of satellites
        for iFreq=1:nFreq
            [labelSatellites,idxBase,idxRover] = intersect(labelSatBase{iFreq},labelSatRover{iFreq});
            SatelliteInformation(iEpoch).obsParametersBase.C{iFreq}  = obsParametersBaseAux.C{iFreq}(idxBase);
            SatelliteInformation(iEpoch).obsParametersBase.L{iFreq}  = obsParametersBaseAux.L{iFreq}(idxBase)*WAVELENGTHS(iFreq);
            SatelliteInformation(iEpoch).posSatBase{iFreq}           = posSatBase{iFreq}(idxBase,:);
            SatelliteInformation(iEpoch).labelSatBase{iFreq}         = labelSatBase{iFreq}(idxBase);
            SatelliteInformation(iEpoch).elevSatBase{iFreq}          = elevSatBase{iFreq}(idxBase);  
                        
            SatelliteInformation(iEpoch).obsParametersRover.C{iFreq} = obsParametersRoverAux.C{iFreq}(idxRover);
            SatelliteInformation(iEpoch).obsParametersRover.L{iFreq} = obsParametersRoverAux.L{iFreq}(idxRover)*WAVELENGTHS(iFreq);
            SatelliteInformation(iEpoch).posSatRover{iFreq}          = posSatRover{iFreq}(idxRover,:);            
            SatelliteInformation(iEpoch).labelSatRover{iFreq}        = labelSatRover{iFreq}(idxRover);  
            SatelliteInformation(iEpoch).elevSatRover{iFreq}         = elevSatRover{iFreq}(idxRover);  
        end              
    end  %End of loading Rover Data    
end %End of Initialization loop



