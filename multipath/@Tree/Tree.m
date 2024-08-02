function new=tree(Params)

%
% Copyright © 2005
% Deutsches Zentrum für Luft- und Raumfahrt e.V.
% German Aerospace Centre
%
% Needed Parameters:
%
% Params.Position           = bottom position of the tree
% Params.Diameter           = diameter of the tree
% Params.ObstacleHeight     = total height of the tree
% Params.MeanDampingFactor  = mean damping factor in dB/m
% Params.TreeTrunkLength    = length of the trunk [m]
% Params.TreeTrunkDiameter  = diameter of the tree trunk [m]
% Params.FresnelRange       = look-up table, or [] for direct calculation
% Params.ScenerySize        = size of the scenery [m] important for display
% Params.GraphicalPlotArea  = visible part of the scenery [m]
% Params.RiceFactor         = Rice Factor of the tree Fading Process
% Params.Bandwidth          = 3 dB Bandwidth of the fading process [1/m]

new=Params;
new.MeanPower=1; % Average Power is 1 of the fading process

% ------------------------------------
%        Internal Data Fields
% ------------------------------------

%Initial Definitions

NrOfRepresentatives=200; %Defines How many carriers are used to aproximate the Spectrum of the fading process

% --------------------- Generating complex Carriers ------------------

% calculating spectrum of the fading process = pdf

NrFreqs=100;

f=linspace(-10*new.Bandwidth,10*new.Bandwidth,NrFreqs);
SpectrumdB=-3/new.Bandwidth*abs(f);
SpectrumLin=10.^(SpectrumdB/10);

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

new=class(new,'Tree');
