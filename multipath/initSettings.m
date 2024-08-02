function settings = initSettings()
%Functions initializes and saves settings. Settings can be edited inside of
%the function, updated from the command line or updated using a dedicated
%GUI - "setSettings".  
%
%All settings are described inside function code.
%
%settings = initSettings()
%
%   Inputs: none
%
%   Outputs:
%       settings     - Receiver settings (a structure). 

%--------------------------------------------------------------------------
%                           SoftGNSS v4.0
% 
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis
% Adapted and maintained by Javier Arribas (jarribas@cttc.es)
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

% CVS record:
% $Id: initSettings.m,v 1.9.2.31 2006/08/18 11:41:57 dpl Exp $

%% Processing settings ====================================================
% Number of milliseconds to be processed used 36000 + any transients (see
% below - in Nav parameters) to ensure nav subframes are provided
settings.msToProcess        = 5000 - 2;        %[ms]

% Number of channels to be used for signal processing
settings.numberOfChannels   = 8;

% Move the starting point of processing. Can be used to start the signal
% processing at any point in the data record (e.g. for long records). fseek
% function is used to move the file read point, therefore advance is byte
% based only. 
settings.skipMs = 1;
settings.samplingFreq       = 20.46e6;
settings.skipNumberOfBytes     = settings.skipMs * settings.samplingFreq * 1e-3;

%% Raw signal file name and other parameter ===============================
% This is a "default" name of the data file (signal record) to be used in
% the post-processing mode
settings.fileName            = '';%'load_data/RxSignal_8sat_20_46MHz_multipath_error_0.dat';%'load_data/RxSignal_6sat_10MHz_multipath.dat';
% settings.fileName            = '/newrawdata.dat';
% settings.fileName            = '/newdata_robustness.dat';
%settings.fileName            = 'signal.dat';
% Data type used to store one sample
%settings.dataType           = 'schar'; % for SiGe v2
%settings.dataType           = 'short'; %for DBF Platform or
% settings.dataType           = 'float'; %For USRP file_sink complex with 32 bits per sample
settings.dataType           = 'float';
% File Types
%1 - 8 bit real samples S0,S1,S2,...
%2 - 8 bit I/Q samples I0,Q0,I1,Q1,I2,Q2,...(SiGe v2)   
%4 - 16 bit I/Q samples I0,Q0,I1,Q1,I2,Q2,...   
%8 - 32 bit I/Q samples I0,Q0,I1,Q1,I2,Q2,... (GNU Radio gr_complex)
settings.fileType           = 8;
% Intermediate, sampling and code frequencies
settings.IF                 = 0;%4.024e6 ;%[Hz]
settings.codeFreqBasis      = 1.023e6;      
settings.GPS_L1             = 1575.42e6; %[Hz]

% Define number of chips in a code period
settings.codeLength         = 1023;

%% Acquisition settings ===================================================
% Skips acquisition in the script postProcessing.m if set to 1
settings.skipAcquisition    = 0;
% List of satellites to look for. Some satellites can be excluded to speed
% up acquisition

%for EGNOS and WAAS satellites (true_PRN = PRN + 87)
%(True PRN are 120 124 126)
settings.acqSatelliteList   = 1:32;         %[PRN numbers] 

% Band around IF to search for satellite signal. Depends on max Doppler
settings.acqSearchBand      = 20;           %[kHz]
% Threshold for the signal presence decision rule
settings.acqThreshold       = 2.5;

%turns on/off the multi-dwell non-coherent integration at acquisition
settings.highsensitivity_acquisition=0;
settings.enable_acquisition_debug_plots=0;
settings.non_coherent_ms=10; % if the highsensitivity_acquisition is enabled (max 10 ms)


%% Tracking loops settings ================================================
% Code tracking loop parameters
settings.dllDampingRatio         = 0.7;
settings.dllNoiseBandwidth       = 2.0;       %[Hz]
settings.dllCorrelatorSpacing    = 0.5;     %[chips]

% Carrier tracking loop parameters
settings.pllDampingRatio         = 0.7;
settings.pllNoiseBandwidth       = 20; %10;%45;      %[Hz]

%% Navigation solution settings ===========================================

% Period for calculating pseudoranges and position
settings.navSolPeriod       = 10;          %[ms]

% Elevation mask to exclude signals from satellites at low elevation
settings.elevationMask      = 0;           %[degrees 0 - 90]
% Enable/dissable use of tropospheric correction
settings.useTropCorr        = 0;            % 0 - Off
                                            % 1 - On
settings.useIonoCorr        = 0;

% True position of the antenna in UTM system (if known). Otherwise enter
% all NaN's and mean position will be used as a reference .
% utmZone = findUtmZone(30.286502, 120.032669);
% [X, Y, Z] = geo2cart(dms2mat(deg2dms(30.286502)), dms2mat(deg2dms(120.032669)), 100, 5); %WGS84
% [E, N, U] = cart2utm(X, Y, Z, utmZone);
% settings.truePosition.E     = E;
% settings.truePosition.N     = N;
% settings.truePosition.U     = U;
X = 4777973.177; Y = 176346.307; Z = 4207663.62;
settings.truePosition.ECEF = [X, Y, Z];%estPos;%
LLA = ecef2lla([X, Y, Z]);
settings.truePosition.LLA = LLA; 
settings.utmZone = findUtmZone(LLA(1), LLA(2));
[E, N, U] = cart2utm(X, Y, Z, settings.utmZone);
settings.truePosition.E     = E;
settings.truePosition.N     = N;
settings.truePosition.U     = U;
%disbled (use the mean position)
% settings.truePosition.E     = nan;
% settings.truePosition.N     = nan;
% settings.truePosition.U     = nan;

%% Plot settings ==========================================================
% Enable/disable plotting of the tracking results for each channel
settings.plotTracking       = 1;            % 0 - Off
                                            % 1 - On

%% Constants ==============================================================

settings.c                  = 299792458;    % The speed of light, [m/s]
settings.startOffset        = 67.727;       %[ms] Initial sign. travel time
