

pth = cd;
addpath ([pth,'\Knife Edge']);
addpath ([pth,'\Parameter Files']);

Parameters.SampFreq = 4e3;                 % Hz
Parameters.MaximumSpeed = 7;               % km/h
Parameters.SatElevation = 30;              % Deg
Parameters.SatAzimut = -45;                % Deg (North == 0, East == 90, South == 180, West == 270)
Parameters.NumberOfSteps = 10;

FontSize=12;

% ---- General Parameters ----

ChannelParams.CarrierFreq = 1.57542e9;     % Hz
ChannelParams.SampFreq = Parameters.SampFreq;
ChannelParams.EnableDisplay = 0;           % 3D visualization is not available in the free version
ChannelParams.EnableCIRDisplay = 1;        % enables CIR display

% ---- Mode Parameters ----

ChannelParams.UserType = 'Pedestrian';
ChannelParams.Surrounding = 'Suburban';
ChannelParams.AntennaHeight = 1;           % m Height of the Antenna
ChannelParams.MinimalPowerdB = -40;        % Echos below this Limit are not initialised

% ---- UserParameters ---

ChannelParams.DistanceFromRoadMiddle = -5.5;% negative: continental (right), positive: England (left)

% ---- Graphics Parameters ---           

ChannelParams.GraphicalPlotArea = 50;      % 
ChannelParams.ViewVector = [-60,20];       % 3D visualization is not available in the free version
ChannelParams.RoadWidth = 8;               %

% --- Building Params ---

ChannelParams.BuildingRow1 = 1;            % logigal to switch Building Row right(heading 0 deg) on
ChannelParams.BuildingRow2 = 1;            % logigal to switch Building Row left (heading 0 deg) on
ChannelParams.BuildingRow1YPosition = -7;  % m
ChannelParams.BuildingRow2YPosition = 7;   % m

ChannelParams.HouseWidthMean = 16;         % m
ChannelParams.HouseWidthSigma = 11;        % m
ChannelParams.HouseWidthMin = 5;           % m
ChannelParams.HouseHeightMin = 3;          % m
ChannelParams.HouseHeightMax = 15;         % m
ChannelParams.HouseHeightMean = 10;        % m
ChannelParams.HouseHeightSigma = 3.6;      % m
ChannelParams.GapWidthMean = 15;           % m
ChannelParams.GapWidthSigma = 22;          % m
ChannelParams.GapWidthMin = 2;             % m
ChannelParams.BuildingGapLikelihood = 0.28;% lin Value

% --- Tree Params ---

ChannelParams.TreeHeight = 7;              % m
ChannelParams.TreeDiameter = 4;            % m
ChannelParams.TreeTrunkLength = 2;         % m
ChannelParams.TreeTrunkDiameter = .2;      % m

ChannelParams.TreeAttenuation = 1.1;       % dB/m

ChannelParams.TreeRow1Use = 1;             % logical switches tree row 1 on
ChannelParams.TreeRow2Use = 1;             % logical switches tree row 2 on

ChannelParams.TreeRow1YPosition = -5;      % m
ChannelParams.TreeRow2YPosition = 5;       % m

ChannelParams.TreeRow1YSigma = 0.5;        % m
ChannelParams.TreeRow2YSigma = 0.5;        % m

ChannelParams.TreeRow1MeanDistance = 40;   % m
ChannelParams.TreeRow2MeanDistance = 20;   % m

ChannelParams.TreeRow1DistanceSigma = 20;  % m
ChannelParams.TreeRow2DistanceSigma = 20;  % m

% --- Pole Params ---

ChannelParams.PoleHeight = 9;              % m
ChannelParams.PoleDiameter = .2;           % m

ChannelParams.PoleRow1Use = 1;             % logical switches Pole row 1 on
ChannelParams.PoleRow2Use = 0;             % logical switches Pole row 2 on

ChannelParams.PoleRow1YPosition = 0;       % m
ChannelParams.PoleRow2YPosition = -5;      % m

ChannelParams.PoleRow1YSigma = 0.5;        % m
ChannelParams.PoleRow2YSigma = 0.5;        % m

ChannelParams.PoleRow1MeanDistance = 40;   % m
ChannelParams.PoleRow2MeanDistance = 40;   % m

ChannelParams.PoleRow1DistanceSigma = 5;   % m
ChannelParams.PoleRow2DistanceSigma = 5;   % m

% ------------ Initial Settings -------------
% - Anything Below here must not be changed -
% -------------------------------------------

Co=2.99e8; % Speed of Light

MaximumPossibleSpeed=Co*Parameters.SampFreq/ChannelParams.CarrierFreq/2; % To fulfill the sampling Theorem
SamplingTime=1/Parameters.SampFreq;

% --- Initialising the channel object ---


TheChannelObject=LandMobileMultipathChannel(ChannelParams);

TimeVec=0;
ComplexOutputVec=[];

% --- Specify power and delay bins for output statistics ---

pwrvec = [0:-1:-30];            % power bins in dB
dlyvec = [0:10e-9:500e-9];      % delay bins in s

PowerDelayProfile(1:length(pwrvec),1:length(dlyvec)) = 0;   % allocate memory
pwrstp = (pwrvec(end)-pwrvec(1))/(length(pwrvec)-1);        % get step size
dlystp = (dlyvec(end)-dlyvec(1))/(length(dlyvec)-1);        % get step size

% --- start simulation ---


