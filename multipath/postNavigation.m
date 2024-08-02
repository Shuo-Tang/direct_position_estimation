function navSolutions = postNavigation(trackResults, settings, satPos, LOSdelay)
%Function calculates navigation solutions for the receiver (pseudoranges,
%positions). At the end it converts coordinates from the WGS84 system to
%the UTM, geocentric or any additional coordinate system.
%
%[navSolutions, eph] = postNavigation(trackResults, settings)
%
%   Inputs:
%       trackResults    - results from the tracking function (structure
%                       array).
%       settings        - receiver settings.
%   Outputs:
%       navSolutions    - contains measured pseudoranges, receiver
%                       clock error, receiver coordinates in several
%                       coordinate systems (at least ECEF and UTM).
%       eph             - received ephemerides of all SV (structure array).

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis with help from Kristin Larson
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
%$Id: postNavigation.m,v 1.1.2.22 2006/08/09 17:20:11 dpl Exp $

nSat = settings.numberOfChannels;
activeChnList = (1: 1: nSat);
subFrameStart = ones(1, nSat);


%% Initialization =========================================================
navSolutions = [];
navSolutions.subFrameStart = subFrameStart;
transmitTime = 0;
%##########################################################################
%#   Do the satellite and receiver position calculations                  #
%##########################################################################

%% Initialization of current measurement ==================================
for currMeasNr = 1:fix((settings.msToProcess - max(subFrameStart)) / ...
                                                     settings.navSolPeriod)
    %all channels enabled!
    %activeChnList=1:1: settings.numberOfChannels;
    
    % Save list of satellites used for position calculation
    navSolutions.channel.PRN(activeChnList, currMeasNr) = ...
                                        [trackResults(activeChnList).PRN]; 

    % These two lines help the skyPlot function. The satellites excluded
    % do to elevation mask will not "jump" to possition (0,0) in the sky
    % plot.
    navSolutions.channel.el(:, currMeasNr) = ...
                                         NaN(settings.numberOfChannels, 1);
    navSolutions.channel.az(:, currMeasNr) = ...
                                         NaN(settings.numberOfChannels, 1);

%% Find pseudoranges ======================================================
     navSolutions.channel.rawP(:, currMeasNr) = calculatePseudoranges(...
             trackResults, ...
             subFrameStart + settings.navSolPeriod * (currMeasNr-1), ...
             activeChnList, settings, LOSdelay);

%% Find satellites positions and clocks corrections =======================
    satPositions = satPos(:,:,(currMeasNr-1)*settings.navSolPeriod + settings.skipMs + 1)';%currMeasNr + settings.skipMs
    % satPositions = satPos(:,:,1)';
    satClkCorr = zeros(1, nSat);
  
   
   navSolutions.transmitTime(currMeasNr)=transmitTime;
   navSolutions.satPositions(currMeasNr,:,:)=satPositions;
   %%%%%%%%%%%%%

%% TODO: (Javi) write RINEX files here
                                    
%% Find receiver position =================================================

    % 3D receiver position can be found only if signals from more than 3
    % satellites are available  
    if length(activeChnList) > 3

        %=== Calculate receiver position ==================================
        %****ORIGINAL ***** 
        [xyzdt, ...
         navSolutions.channel.el(activeChnList, currMeasNr), ...
         navSolutions.channel.az(activeChnList, currMeasNr), ...
         navSolutions.DOP(:, currMeasNr), ...
         navSolutions.Q_DOP(:,:,currMeasNr)] = ...
            leastSquarePos(satPositions, ...
                           navSolutions.channel.rawP(activeChnList, currMeasNr)' + satClkCorr * settings.c, ...
                           settings);
                               

        %--- Save results -------------------------------------------------

        %----------
        navSolutions.X(currMeasNr)  = xyzdt(1);
        navSolutions.Y(currMeasNr)  = xyzdt(2);
        navSolutions.Z(currMeasNr)  = xyzdt(3);
        navSolutions.dt(currMeasNr) = xyzdt(4);

        %=== Correct pseudorange measurements for clocks errors ===========
         navSolutions.channel.correctedP(activeChnList, currMeasNr) = ...
                 navSolutions.channel.rawP(activeChnList, currMeasNr) + ...
                 satClkCorr' * settings.c + navSolutions.dt(currMeasNr);            

%% Coordinate conversion ==================================================

        %=== Convert to geodetic coordinates ==============================
        [navSolutions.latitude(currMeasNr), ...
         navSolutions.longitude(currMeasNr), ...
         navSolutions.height(currMeasNr)] = cart2geo(...
                                            navSolutions.X(currMeasNr), ...
                                            navSolutions.Y(currMeasNr), ...
                                            navSolutions.Z(currMeasNr), ...
                                            5);

        %=== Convert to UTM coordinate system =============================
        navSolutions.utmZone = findUtmZone(navSolutions.latitude(currMeasNr), ...
                                           navSolutions.longitude(currMeasNr)); 
        
        [navSolutions.E(currMeasNr), ...
         navSolutions.N(currMeasNr), ...
         navSolutions.U(currMeasNr)] = cart2utm(xyzdt(1), xyzdt(2), ...
                                                xyzdt(3), ...
                                                navSolutions.utmZone);
        
    else % if size(activeChnList, 2) > 3 
        %--- There are not enough satellites to find 3D position ----------
        disp(['   Measurement No. ', num2str(currMeasNr), ...
                       ': Not enough information for position solution.']);

        %--- Set the missing solutions to NaN. These results will be
        %excluded automatically in all plots. For DOP it is easier to use
        %zeros. NaN values might need to be excluded from results in some
        %of further processing to obtain correct results.
        navSolutions.X(currMeasNr)           = NaN;
        navSolutions.Y(currMeasNr)           = NaN;
        navSolutions.Z(currMeasNr)           = NaN;
        navSolutions.dt(currMeasNr)          = NaN;
        navSolutions.DOP(:, currMeasNr)      = zeros(5, 1);
        navSolutions.latitude(currMeasNr)    = NaN;
        navSolutions.longitude(currMeasNr)   = NaN;
        navSolutions.height(currMeasNr)      = NaN;
        navSolutions.E(currMeasNr)           = NaN;
        navSolutions.N(currMeasNr)           = NaN;
        navSolutions.U(currMeasNr)           = NaN;

        navSolutions.channel.az(activeChnList, currMeasNr) = ...
                                             NaN(1, length(activeChnList));
        navSolutions.channel.el(activeChnList, currMeasNr) = ...
                                             NaN(1, length(activeChnList));

        % TODO: Know issue. Satellite positions are not updated if the
        % satellites are excluded do to elevation mask. Therefore rasing
        % satellites will be not included even if they will be above
        % elevation mask at some point. This would be a good place to
        % update positions of the excluded satellites.

    end % if size(activeChnList, 2) > 3

    %=== Update the transmit time ("measurement time") ====================
    transmitTime = transmitTime + settings.navSolPeriod / 1000;

end %for currMeasNr...
