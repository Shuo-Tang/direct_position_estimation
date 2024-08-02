function new=LandMobileMultipathChannel(Params)

%
% Copyright © 2005
% Deutsches Zentrum für Luft- und Raumfahrt e.V.
% German Aerospace Centre
%
% This is the Master Module for the Land mobile multipath channel.
%
% call: new=LandMobileMultipathChannel(ParameterStruct)
%

ClassName='LandMobileMultipathChannel';
MethodName='LandMobileMultipathChannel';

% Tree Fading Process

TreeBandwidth       = 0.437; %[1/m] = 3 dB Bandwidth of the tree fading spectrum
TreeRiceFactor      = 30; 

% create a look-up table for Fresnel range
FresnelRangeTable=FresnelRange(0,20,0.01);

% Security Checks for the Parameters
%

% ------------ MasterSwitch --------------

if (isequal(Params.Surrounding,'Urban') | isequal(Params.Surrounding,'Suburban')) & (isequal(Params.UserType,'Car') | isequal(Params.UserType,'Pedestrian'))

    Params.MaximumNumberOfReflectors=80;

else
    Block_error(ClassName,MethodName,'Sorry only Urban and Suburban Car or Pedestrian channels are implemented at the moment')
end%if

if Params.SampFreq<0
    Block_error(ClassName,MethodName,'Parameter Sampling Frequency must be > 0')
end%if

if ~any(Params.EnableDisplay==[0,1])
    Block_error(ClassName,MethodName,'Parameter EnableDisplay must be 0 or 1')
end%if



% ----- taking over the parameters -----

new=Params;
new.SatVector=[];

% ------------ Declaring Member Fields ------------

new.CurrentRxPosition=[0 new.DistanceFromRoadMiddle new.AntennaHeight];
new.Time=0;
new.DisplayAltitude=50;     %m Defines the display Altitude
new.MaskingAngle=5;         %Deg Below this Angle no Obstacles are generated
new.LastLOSPhase=0;         % Store laste phase of the line of sight

% ---------- internalCalculations -----------
Co=2.99e8;
new.WaveLength=Co/new.CarrierFreq;


% calculating the distance of the largest obstacle not to interfere with the masking angle
MaximumObstacleHeight = max([new.HouseHeightMax,new.TreeHeight,new.PoleHeight]);
new.ScenerySize=(MaximumObstacleHeight+10*new.WaveLength)/tan(new.MaskingAngle/180*pi);

disp(['Scenery Size ',num2str(new.ScenerySize)])
% Initialising the sub modules

% ----- Reflectors -----

ReflectorListParams.MaximumNumberOfReflectors=Params.MaximumNumberOfReflectors;
ReflectorListParams.EnableDisplay=Params.EnableDisplay;
ReflectorListParams.CarrierFreq=Params.CarrierFreq;

new.TheReflectorList=ListOfReflectors(ReflectorListParams);


% ---- Tree Row 1 ----

if isequal(new.TreeRow1Use,1)
    Tree1Parameters.ScenerySize                 = new.ScenerySize;
    Tree1Parameters.Content                     ='Tree';
    Tree1Parameters.ObstacleHeight              = new.TreeHeight;
    Tree1Parameters.ObstacleDiameter            = new.TreeDiameter;
    Tree1Parameters.ObstacleYPosition           = new.TreeRow1YPosition;
    Tree1Parameters.ObstacleYSigma              = new.TreeRow1YSigma;
    Tree1Parameters.ObstacleMeanDistance        = new.TreeRow1MeanDistance;
    Tree1Parameters.ObstacleDistanceSigma       = new.TreeRow1DistanceSigma;
    Tree1Parameters.TreeAttenuation             = new.TreeAttenuation;
    Tree1Parameters.TreeTrunkLength             = new.TreeTrunkLength;
    Tree1Parameters.TreeTrunkDiameter           = new.TreeTrunkDiameter;
    Tree1Parameters.FresnelRange                = FresnelRangeTable;
    Tree1Parameters.GraphicalPlotArea           = new.GraphicalPlotArea;
    Tree1Parameters.TreeBandwidth               = TreeBandwidth;
    Tree1Parameters.TreeRiceFactor              = TreeRiceFactor;
    new.TreeRow1=ListOfObstacles(Tree1Parameters);
    new.TreeRow1=CreateNewKillOld(new.TreeRow1,new.CurrentRxPosition);
end%if

% ---- Tree Row 2 ----

if isequal(new.TreeRow2Use,1)
    Tree2Parameters.ScenerySize                 = new.ScenerySize;
    Tree2Parameters.Content                     ='Tree';
    Tree2Parameters.ObstacleHeight              = new.TreeHeight;
    Tree2Parameters.ObstacleDiameter            = new.TreeDiameter;
    Tree2Parameters.ObstacleYPosition           = new.TreeRow2YPosition;
    Tree2Parameters.ObstacleYSigma              = new.TreeRow2YSigma;
    Tree2Parameters.ObstacleMeanDistance        = new.TreeRow2MeanDistance;
    Tree2Parameters.ObstacleDistanceSigma       = new.TreeRow2DistanceSigma;
    Tree2Parameters.TreeAttenuation             = new.TreeAttenuation;
    Tree2Parameters.TreeTrunkLength             = new.TreeTrunkLength;
    Tree2Parameters.TreeTrunkDiameter           = new.TreeTrunkDiameter;
    Tree2Parameters.FresnelRange                = FresnelRangeTable;
    Tree2Parameters.GraphicalPlotArea           = new.GraphicalPlotArea;
    Tree2Parameters.TreeBandwidth               = TreeBandwidth;
    Tree2Parameters.TreeRiceFactor              = TreeRiceFactor;


    new.TreeRow2=ListOfObstacles(Tree2Parameters);
    new.TreeRow2=CreateNewKillOld(new.TreeRow2,new.CurrentRxPosition);

end%if

% ---- Pole Row 1 ----

if isequal(new.PoleRow1Use,1)
    Pole1Parameters.ScenerySize                 = new.ScenerySize;
    Pole1Parameters.Content                     ='Pole';
    Pole1Parameters.ObstacleHeight              = new.PoleHeight;
    Pole1Parameters.ObstacleDiameter            = new.PoleDiameter;
    Pole1Parameters.ObstacleYPosition           = new.PoleRow1YPosition;
    Pole1Parameters.ObstacleYSigma              = new.PoleRow1YSigma;
    Pole1Parameters.ObstacleMeanDistance        = new.PoleRow1MeanDistance;
    Pole1Parameters.ObstacleDistanceSigma       = new.PoleRow1DistanceSigma;
    Pole1Parameters.TreeAttenuation             = new.TreeAttenuation;
    Pole1Parameters.TreeTrunkLength             = new.TreeTrunkLength;
    Pole1Parameters.TreeTrunkDiameter           = new.TreeTrunkDiameter;
    Pole1Parameters.FresnelRange                = FresnelRangeTable;
    Pole1Parameters.GraphicalPlotArea           = new.GraphicalPlotArea;
    Pole1Parameters.TreeBandwidth               = TreeBandwidth;
    Pole1Parameters.TreeRiceFactor              = TreeRiceFactor;



    new.PoleRow1=ListOfObstacles(Pole1Parameters);
    new.PoleRow1=CreateNewKillOld(new.PoleRow1,new.CurrentRxPosition);
end%if

% ---- Pole Row 2 ----

if isequal(new.PoleRow2Use,1)
    Pole2Parameters.ScenerySize                 = new.ScenerySize;
    Pole2Parameters.Content                     ='Pole';
    Pole2Parameters.ObstacleHeight              = new.PoleHeight;
    Pole2Parameters.ObstacleDiameter            = new.PoleDiameter;
    Pole2Parameters.ObstacleYPosition           = new.PoleRow2YPosition;
    Pole2Parameters.ObstacleYSigma              = new.PoleRow2YSigma;
    Pole2Parameters.ObstacleMeanDistance        = new.PoleRow2MeanDistance;
    Pole2Parameters.ObstacleDistanceSigma       = new.PoleRow2DistanceSigma;
    Pole2Parameters.TreeAttenuation             = new.TreeAttenuation;
    Pole2Parameters.TreeTrunkLength             = new.TreeTrunkLength;
    Pole2Parameters.TreeTrunkDiameter           = new.TreeTrunkDiameter;
    Pole2Parameters.FresnelRange                = FresnelRangeTable;
    Pole2Parameters.GraphicalPlotArea           = new.GraphicalPlotArea;
    Pole2Parameters.TreeBandwidth               = TreeBandwidth;
    Pole2Parameters.TreeRiceFactor              = TreeRiceFactor;



    new.PoleRow2=ListOfObstacles(Pole2Parameters);
    new.PoleRow2=CreateNewKillOld(new.PoleRow2,new.CurrentRxPosition);
end%if

if new.BuildingRow1
    %  Housefront 1

    HouseFront1Params.ScenerySize       = new.ScenerySize;
    HouseFront1Params.HouseWidthMean    = new.HouseWidthMean;
    HouseFront1Params.HouseWidthSigma   = new.HouseWidthSigma;
    HouseFront1Params.HouseHeightMean   = new.HouseHeightMean;
    HouseFront1Params.HouseHeightSigma  = new.HouseHeightSigma;
    HouseFront1Params.HouseHeightMax    = new.HouseHeightMax;
    HouseFront1Params.GraphicalPlotArea = new.GraphicalPlotArea;
    HouseFront1Params.GapLikelihood     = new.BuildingGapLikelihood;
    HouseFront1Params.HouseWidthMin     = new.HouseWidthMin;
    HouseFront1Params.HouseWidthMin     = new.HouseWidthMin;
    HouseFront1Params.HouseHeightMin    = new.HouseHeightMin;
    HouseFront1Params.HouseHeightMax    = new.HouseHeightMax;
    HouseFront1Params.GapWidthMean      = new.GapWidthMean;
    HouseFront1Params.GapWidthSigma     = new.GapWidthSigma;
    HouseFront1Params.GapWidthMin       = new.GapWidthMin;

    % Initialize HouseFront 1

    HouseFront1Params.YPosition         = new.BuildingRow1YPosition;
    new.HouseFront1=Housefronts(HouseFront1Params);
end%if

if new.BuildingRow2
    %  Housefront 1

    HouseFront2Params.ScenerySize       = new.ScenerySize;
    HouseFront2Params.HouseWidthMean    = new.HouseWidthMean;
    HouseFront2Params.HouseWidthSigma   = new.HouseWidthSigma;
    HouseFront2Params.HouseHeightMean   = new.HouseHeightMean;
    HouseFront2Params.HouseHeightSigma  = new.HouseHeightSigma;
    HouseFront2Params.HouseHeightMax    = new.HouseHeightMax;
    HouseFront2Params.GraphicalPlotArea = new.GraphicalPlotArea;
    HouseFront2Params.GapLikelihood     = new.BuildingGapLikelihood;
    HouseFront2Params.HouseWidthMin     = new.HouseWidthMin;
    HouseFront2Params.HouseWidthMin     = new.HouseWidthMin;
    HouseFront2Params.HouseHeightMin    = new.HouseHeightMin;
    HouseFront2Params.HouseHeightMax    = new.HouseHeightMax;
    HouseFront2Params.GapWidthMean      = new.GapWidthMean;
    HouseFront2Params.GapWidthSigma     = new.GapWidthSigma;
    HouseFront2Params.GapWidthMin       = new.GapWidthMin;

    HouseFront2Params.YPosition         = new.BuildingRow2YPosition;

    % Initialize HouseFront 2

    HouseFront2Params.YPosition         = new.BuildingRow2YPosition;
    new.HouseFront2=Housefronts(HouseFront2Params);
end%if

% ---- Initialise Echo Number Generator  &  Echo Generator ---

switch Params.UserType
    case 'Car'
        switch Params.Surrounding
            case 'Urban'
                new.NumberGenerator=EchoNumberGenerator('EchoNumberParUrbanCar.mat');
                new.RandomEchoGenerator=EchoGenerator('EchoParUrbanCar.mat');
            case 'Suburban'
                new.NumberGenerator=EchoNumberGenerator('EchoNumberParSuburbanCar.mat');
                new.RandomEchoGenerator=EchoGenerator('EchoParSuburbanCar.mat');
            otherwise
                error('Echo Generator data not present for this channel type')
        end
    case 'Pedestrian'
        switch Params.Surrounding
            case 'Urban'
                new.NumberGenerator=EchoNumberGenerator('EchoNumberParUrbanPedestrian.mat');
                new.RandomEchoGenerator=EchoGenerator('EchoParUrbanPedestrian.mat');
            case 'Suburban'
                new.NumberGenerator=EchoNumberGenerator('EchoNumberParSuburbanPedestrian.mat');
                new.RandomEchoGenerator=EchoGenerator('EchoParSuburbanPedestrian.mat');
            otherwise
                error('Echo Generator data not present for this channel type')
        end
    otherwise
        Block_error(ClassName,MethodName,'Sorry only Urban and Suburban Car and Pedestrian channels are implemented at the moment')
end

new=class(new,'LandMobileMultipathChannel');



















