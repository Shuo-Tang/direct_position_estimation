function new=Pole(Params)

%
% Copyright © 2005
% Deutsches Zentrum für Luft- und Raumfahrt e.V.
% German Aerospace Centre
%
% Needed Parameters:
%
% Params.Position           = bottom position of the pole
% Params.Diameter           = diameter of the pole
% Params.ObstacleHeight     = height of the pole
% Params.FresnelRange       = look-up table, or [] for direct calculation
% Params.ScenerySize        = size of the scenery [m] important for display
% Params.GraphicalPlotArea  = visible part of the scenery [m]

new=Params;

new=class(new,'Pole');
