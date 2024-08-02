function obout=Reflector(varargin)

%
% Copyright © 2005
% Deutsches Zentrum für Luft- und Raumfahrt e.V.
% German Aerospace Centre
%
% Constructor of class Reflector
%
% function call:
%
% Reflector(ParameterStruct)
%
% ------------------------------------
%        Parameter data Fields
% ------------------------------------
%
%         ParameterStruct.Position                  % [x y z] Coordinates of the reflector [m]
%         ParameterStruct.Bandwidth                 % 3dB Bandwidth of the reflection in Hz (two side band)
%         ParameterStruct.MeanPower                 % Mean Power of the reflector [dB]
%         ParameterStruct.Phase                     % Initial Phase of the reflector [rad]
%         ParameterStruct.RiceFactor                % Defines the rice factor of the reflector
%         ParameterStruct.EnableDisplay             % True if display is enabled
%         ParameterStruct.CarrierFreq               % CarrierFrequency in Hz
%         ParameterStruct.SpeedRatio                % gives the ratio
%         between the receiver speed and the speed of the Reflector. Usually 0 or 1

ClassName  = 'Reflector';
MethodName = 'Reflector';

switch nargin
    case 0
        Block_error(ClassName,MethodName,'Contructor called without arguments')
    case 1


        
        % ------------------------------------
        % making the Parameters become members:
        % ------------------------------------
        
        new=varargin{1};

        %calculating linear power from dB
        
         new.MeanPower=10^(new.MeanPower/10);
         
        % ------------------------------------
        %        Internal Data Fields
        % ------------------------------------

        %Initial Definitions

        NrOfRepresentatives=200; %Defines How many carriers are used to aproximate the Spectrum of the fading process
        Co=2.99e8;
        new.Wavelength=Co/new.CarrierFreq;
        % --------- filling the fields --------------

        if new.EnableDisplay
            new.Display.ScalingFactor=7;        % defines the ratio between Power and ball size in meters
            new.Display.NrOfSpherePoints=30;    % defines how many points are used to plot the sphere
            new.Display.NrOfColorSteps=128;     % defines how many colorsteps are used for display
            new.Display.MinPower=-40;           % [dB] defines the minimal Power used for the color coding
            new.Display.MaxPower=-10;           % [dB] defines the maximal Power used for the color coding
            new.Display.MinBallSize=.1;         % Defines the minimun´m Ball size when the Power is <= MinPower
            % --------- filling the fields --------------
            
            new.Display.Colormap=jet(new.Display.NrOfColorSteps);
            [new.Display.UnitSphereX,new.Display.UnitSphereY,new.Display.UnitSphereZ] = sphere(new.Display.NrOfSpherePoints);
            new.Display.PowerLookupVector=linspace(0,1,new.Display.NrOfColorSteps);
        else
            new.Display.ScalingFactor=[];
            new.Display.NrOfSpherePoints=[];
            new.Display.NrOfColorSteps=[];
            new.Display.MinPower=[];
            new.Display.MaxPower=[];
            new.Display.MinBallSize=[];
        end% if


         new.ActualPowerdB=[];
        
        
        % --------------------- Generating complex Carriers ------------------

        % calculating spectrum of the fading process = pdf
        
        NrFreqs=100;
        SigmaOfSpectrum = new.Bandwidth/2/sqrt(log(2)*2); % 3dB Bandwidth therefore the scaling factor sqrt(log(2)*2), Two processes are generated therefore the scaling factor 1/2

        f=linspace(-3*SigmaOfSpectrum,3*SigmaOfSpectrum,NrFreqs);
        SpectrumLin= 1/sqrt(2*pi*SigmaOfSpectrum^2)*exp(-f.^2/2/SigmaOfSpectrum^2); %Gaussian Function

        % calculating cumulative distribution function (cdf)
        
        
        cdf=cumsum(SpectrumLin);
        cdf=cdf/cdf(end);

        % Randomly variable [0...1]
        
        UniRandom=rand(1,NrOfRepresentatives);
        
        % Calculating the Frequencies by inverse table lookup
        
        new.RandomFrequencies=interp1(cdf,f,UniRandom,'spline');
       

        % ----- calculation of the Phases -----
        
        new.RandomPhases=2*pi*rand(1,NrOfRepresentatives);

        % ----- calculation of the parameters for the rician distribution -----

        % RiceFactor new.RiceFactor      
        
        L=sqrt(2*new.RiceFactor);

        % Power of a Std Rice Process

        P=2+L^2;

        new.MyRiceanConst=L*exp(j*rand(1,1)*2*pi);

        new.ProcessScalingFactor=1/sqrt(P)*sqrt(new.MeanPower);

        
        
        
        
        
        
    otherwise
            Block_error(ClassName,MethodName,'Contructor called with the wrong number of arguments. 1 Paramter is required - the parameter struct')
end

obout=class(new,'Reflector');

