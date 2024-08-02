function this = CreateNewKillOld(this,CurrentPosition)

%
% Copyright © 2005
% Deutsches Zentrum für Luft- und Raumfahrt e.V.
% German Aerospace Centre
%
% determines Obstacles which had left the Scenery and creates new ones if
% needed
%

% create new ones until the cenery end is reached

while this.ObstacleCoordinate{end}* [1 0 0]' - CurrentPosition(1) < this.ScenerySize
    Distance=this.ObstacleMeanDistance+this.ObstacleDistanceSigma*randn(1,1);
    YVariation=this.ObstacleYSigma*randn(1,1);
    YVariation=max(-this.ScenerySize,YVariation);
    YVariation=min(this.ScenerySize,YVariation);

    LastXCoordinate=this.ObstacleCoordinate{end}*[1,0,0]';

    this.ObstacleCoordinate{end+1}=[LastXCoordinate+Distance,this.ObstacleYPosition+YVariation,this.ZCoodinate];
    switch this.Content
        case'Tree'
            ObstacleParams=[];
            ObstacleParams.Position             = this.ObstacleCoordinate{end};
            ObstacleParams.Diameter             = this.ObstacleDiameter;
            ObstacleParams.MeanDampingFactor    = this.TreeAttenuation;
            ObstacleParams.Bandwidth            = this.TreeBandwidth;
            ObstacleParams.RiceFactor           = this.TreeRiceFactor;
            ObstacleParams.ObstacleHeight       = this.ObstacleHeight;
            ObstacleParams.FresnelRange         = this.FresnelRange;
            ObstacleParams.GraphicalPlotArea    = this.GraphicalPlotArea;
            ObstacleParams.TreeTrunkLength      = this.TreeTrunkLength;
            ObstacleParams.TreeTrunkDiameter    = this.TreeTrunkDiameter;




            this.TheList{end+1}                 = Tree(ObstacleParams);

        case 'Pole'
            ObstacleParams=[];
            ObstacleParams.Position             = this.ObstacleCoordinate{end};
            ObstacleParams.Diameter             = this.ObstacleDiameter;
            ObstacleParams.ObstacleHeight       = this.ObstacleHeight;
            ObstacleParams.FresnelRange         = this.FresnelRange;
            ObstacleParams.GraphicalPlotArea    = this.GraphicalPlotArea;


            this.TheList{end+1}                 = Pole(ObstacleParams);
        otherwise
            error('Content not implemented')
    end%switch
end %while

% kill those which have left the Scenery

while this.ObstacleCoordinate{1}* [1 0 0]' - CurrentPosition(1) < -this.ScenerySize

    this.TheList=this.TheList(2:end);
    this.ObstacleCoordinate=this.ObstacleCoordinate(2:end);

end%while


