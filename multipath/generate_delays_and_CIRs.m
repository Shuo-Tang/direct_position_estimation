function [LOS_amplitudes, LOS_delays, echo_amplitudes, echo_delays] = ...
    generate_delays_and_CIRs(user_type, scenario, Parameters)
% clear classes
% LandMobileMultipathChannel_Demo_UrbanCar
% Version 3.0
%% ---- General Parameters ----
% === irrelevant parameters (but must include) === % 
ChannelParams.EnableDisplay=0;           % 3D visualization is not available in the free version
ChannelParams.EnableCIRDisplay=1;        % enables CIR display

%% ---- Mode Parameters ----
% === common parameters === %
% channel
ChannelParams.CarrierFreq = 1.57542e9;     % Hz
ChannelParams.SampFreq = Parameters.SampFreq;
ChannelParams.MinimalPowerdB = -40;        % Echos below this Limit are not initialised
% visualization
ChannelParams.GraphicalPlotArea = 50;      % 
ChannelParams.ViewVector = [-60,20];     % 3D visualization is not available in the free version
% building
ChannelParams.BuildingRow1 = 1;            % logigal to switch Building Row right(heading 0 deg) on
ChannelParams.BuildingRow2 = 1;            % logigal to switch Building Row left (heading 0 deg) on

% === classified parameters === %
if strcmp(user_type, 'Car')
    ChannelParams.UserType = 'Car';
    ChannelParams.AntennaHeight = 2;         % m Height of the Antenna
elseif strcmp(user_type, 'Pedestrian')
    ChannelParams.UserType = 'Pedestrian';
    ChannelParams.AntennaHeight = 1;         % m Height of the Antenna
end

if strcmp(scenario, 'Urban')
    ChannelParams.Surrounding = 'Urban';
    % road
    ChannelParams.RoadWidth = 15;
    % building
    ChannelParams.BuildingRow1YPosition = -12; % m
    ChannelParams.BuildingRow2YPosition = 12;  % m
    ChannelParams.HouseWidthMean = 22;         % m
    ChannelParams.HouseWidthSigma = 25;        % m
    ChannelParams.HouseWidthMin = 10;          % m
    ChannelParams.HouseHeightMin = 4;          % m
    ChannelParams.HouseHeightMax  =  50;         % m
    ChannelParams.HouseHeightMean = 16;        % m
    ChannelParams.HouseHeightSigma = 6.4;      % m
    ChannelParams.GapWidthMean = 27;           % m
    ChannelParams.GapWidthSigma = 25;          % m
    ChannelParams.GapWidthMin = 10;            % m
    ChannelParams.BuildingGapLikelihood = 0.18;% lin Value
    % tree
    ChannelParams.TreeHeight = 6;            % m
    ChannelParams.TreeDiameter = 3;          % m
    ChannelParams.TreeTrunkLength = 2;         % m
    ChannelParams.TreeTrunkDiameter = .2;      % m
    ChannelParams.TreeAttenuation = 1.1;     % dB/m
    ChannelParams.TreeRow1Use = 1;             % logical switches tree row 1 on
    ChannelParams.TreeRow2Use = 1;             % logical switches tree row 2 on
    ChannelParams.TreeRow1YPosition = -6;      % m
    ChannelParams.TreeRow2YPosition = 6;       % m
    ChannelParams.TreeRow1YSigma = 0.5;        % m
    ChannelParams.TreeRow2YSigma = 0.5;        % m
    ChannelParams.TreeRow1MeanDistance = 60;   % m
    ChannelParams.TreeRow2MeanDistance = 40;   % m
    ChannelParams.TreeRow1DistanceSigma = 20;  % m
    ChannelParams.TreeRow2DistanceSigma = 20;  % m
    % pole
    ChannelParams.PoleHeight = 10;           % m
    ChannelParams.PoleDiameter = .2;         % m
    ChannelParams.PoleRow1Use = 1;             % logical switches Pole row 1 on
    ChannelParams.PoleRow2Use = 0;             % logical switches Pole row 2 on
    ChannelParams.PoleRow1YPosition = 0;       % m
    ChannelParams.PoleRow2YPosition = 10;      % m
    ChannelParams.PoleRow1YSigma = 1;          % m
    ChannelParams.PoleRow2YSigma = 1;          % m
    ChannelParams.PoleRow1MeanDistance = 25;   % m
    ChannelParams.PoleRow2MeanDistance = 10;   % m
    ChannelParams.PoleRow1DistanceSigma = 10;  % m
    ChannelParams.PoleRow2DistanceSigma = 10;  % m
elseif strcmp(scenario, 'Suburban')
    ChannelParams.Surrounding = 'Suburban';
    % road
    ChannelParams.RoadWidth = 8; 
    % building
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
    % tree
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
    % pole
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

end

% === specific parameters === %
if strcmp(user_type, 'Car') && strcmp(scenario, 'Urban')
    ChannelParams.DistanceFromRoadMiddle = -5; % negative: continental (right), positive: England (left)
    ChannelParams.RoadWidth = 15;
    ChannelParams.BuildingRow1YPosition = -12; % m
    ChannelParams.BuildingRow2YPosition = 12;  % m
    ChannelParams.TreeHeight = 8;            % m
    ChannelParams.TreeDiameter = 5;          % m
    ChannelParams.TreeRow1YPosition = -8;      % m
    ChannelParams.TreeRow2YPosition = 8;       % m
    ChannelParams.TreeRow1YSigma = 2;          % m
    ChannelParams.TreeRow2YSigma = 2;          % m
    ChannelParams.PoleRow1YPosition = 0;       % m
    ChannelParams.PoleRow2YPosition = 10;      % m
    ChannelParams.PoleRow1YSigma = 1;          % m
    ChannelParams.PoleRow2YSigma = 1;          % m
end

if strcmp(user_type, 'Car') && strcmp(scenario, 'Suburban')
    ChannelParams.DistanceFromRoadMiddle = -2; 
end

if strcmp(user_type, 'Pedestrian') && strcmp(scenario, 'Urban')
    ChannelParams.DistanceFromRoadMiddle = -6.5;
    ChannelParams.RoadWidth = 10; 
    ChannelParams.BuildingRow1YPosition = -8;  % m
    ChannelParams.BuildingRow2YPosition = 8;   % m
    ChannelParams.TreeHeight = 6;            % m
    ChannelParams.TreeDiameter = 3;          % m
    ChannelParams.TreeRow1YPosition = -6;      % m
    ChannelParams.TreeRow2YPosition = 6;       % m
    ChannelParams.TreeRow1YSigma = 0.5;        % m
    ChannelParams.TreeRow2YSigma = 0.5;        % m
    ChannelParams.PoleRow1YPosition = -6;      % m
    ChannelParams.PoleRow2YPosition = 0;       % m
    ChannelParams.PoleRow1YSigma = 0.5;        % m
    ChannelParams.PoleRow2YSigma = 0.5;        % m
end

if strcmp(user_type, 'Pedestrian') && strcmp(scenario, 'Suburban')
    ChannelParams.DistanceFromRoadMiddle = -5.5;
end



%% ---- Initial Settings -----
% --- Initialising the channel object ---
disp('Initialising the channel ...')
TheChannelObject=LandMobileMultipathChannel(ChannelParams);

TimeVec=0;
ComplexOutputVec=[];

%% --- start simulation ---
% memory allocation
% LOS_coeffs, LOS_delays, EchoCoeffs, EchoDelays
% start simulation
LOS_amplitudes = cell(Parameters.NumberOfSteps, 1);
LOS_delays = cell(Parameters.NumberOfSteps, 1);
echo_amplitudes = cell(Parameters.NumberOfSteps, 1);
echo_delays = cell(Parameters.NumberOfSteps, 1);

for dhv=1:Parameters.NumberOfSteps
    % fprintf('%d \n', dhv);

    TimeVec(end)=dhv/Parameters.SampFreq;

    % --- "drunken" driver movement example ---
    
    ActualSpeed=Parameters.MaximumSpeed/2/3.6*(1+sin(TimeVec(end)/3));   
    % SpeedVec(dhv)=ActualSpeed;              % m/s
    ActualHeading=20*sin(TimeVec(end)/3);   % Deg (North == 0, East == 90, South == 180, West == 270)
    
    % --- generate CIR ---

    [TheChannelObject,LOS,LOSDelays,ComplexOutputVec,DelayVec,EchoNumberVec,WayVec(dhv),TimeVec(dhv)]=generate(TheChannelObject,ActualSpeed,ActualHeading,Parameters.SatElevation,Parameters.SatAzimut);
    % [ChannelObject, LOSCoeff, LOSDelays,...
    %     EchoCoeff, EchoDelays, EchoNumbers,... 
    %     WayVec, TimeVec]= generate( ChannelObject, ActualSpeed, ActualHeading,...
    %                                 SatElevation, SatAzimuth); 
    LOS_amplitudes{dhv} = LOS;
    LOS_delays{dhv} = LOSDelays;
    echo_amplitudes{dhv} = ComplexOutputVec;
    echo_delays{dhv} = DelayVec;
end
1;
