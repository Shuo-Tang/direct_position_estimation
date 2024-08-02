function new=ListOfObstacles(Params)

%
% Copyright © 2005
% Deutsches Zentrum für Luft- und Raumfahrt e.V.
% German Aerospace Centre
%
% call MyObject=ListOfObstacles(Params)
%
% params.ScenerySize   Defines the Scenery Size
% params.Content       Defines what to generate [char]
% params.ObstacleHeight = 10;%m
% params.ObstacleDiameter = 4;%m
% params.ObstacleYPosition=-10; %m
% params.ObstacleYSigma=1; %m
% params.ObstacleMeanDistance=10;%m
% params.ObstacleDistanceSigma=5;%m
% params.GraphicalPlotArea [m] defines the size of the plot area
% params.FresnelRange = look-up table, or [] for direct calculation
% params.TreeBandwidth [1/m] = 3 dB Bandwidth of the tree fading spectrum
% params.RiceFactor = Rice Factor of the Tree fading process

% Initialising the fields

new=Params;
new.ZCoodinate=0; % Defines the lower Z coordinate of the obstacle

new.TheList={};
new.ObstacleCoordinate{1}=[-new.ScenerySize,new.ObstacleYPosition,new.ZCoodinate];

new=class(new,'ListOfObstacles');


