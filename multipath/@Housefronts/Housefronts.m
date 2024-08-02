function new=Housfronts(Params);

%
% Copyright © 2005
% Deutsches Zentrum für Luft- und Raumfahrt e.V.
% German Aerospace Centre
%
% call MyObject=Housefronts(Params)
%
% Params.ScenerySize             
% Params.YPosition
% Params.HouseWidthMean
% Params.HouseWidthSigma
% Params.HouseWidthMin
% Params.HouseHeightMean
% Params.HouseHeightSigma
% Params.HouseHeightMin
% Params.HouseHeightMax
% Params.GapWidthMean
% Params.GapWidthSigma
% Params.GapWidthMin
% Params.GraphicalPlotArea       
% Params.GapLikelihood


new=Params;

new.XValues=[];        % positions of walls (to the right of each house)
new.HeightValues=[];   % height of houses (i.e., to the left of XValue)
new.IDValuesX=[];      % identifier for each wall
new.IDValuesZ=[];      % identifier for each roof
new.IDCounter=1;

new=class(new,'Housefronts');
