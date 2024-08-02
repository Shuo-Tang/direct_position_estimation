
% LandMobileMultipathChannel_Demo_UrbanCar
% Version 3.0

close all
clear all
clear classes
clc

pth = cd;
addpath ([pth,'\Knife Edge']);
addpath ([pth,'\Parameter Files']);
Parameters.SampFreq=10;        % Hz
Parameters.MaximumSpeed=50;     % km/h
Parameters.SatElevation=30;     % Deg
Parameters.SatAzimut=-45;       % Deg (North == 0, East == 90, South == 180, West == 270)
Parameters.NumberOfSteps=600;

FontSize=12;

% ---- General Parameters ----

ChannelParams.CarrierFreq=1.57542e9;     % Hz
ChannelParams.SampFreq=Parameters.SampFreq;
ChannelParams.EnableDisplay=0;           % 3D visualization is not available in the free version
ChannelParams.EnableCIRDisplay=1;        % enables CIR display

% ---- Mode Parameters ----

ChannelParams.UserType = 'Car'; 
ChannelParams.Surrounding = 'Urban'; 
ChannelParams.AntennaHeight = 2;         % m Height of the Antenna1
ChannelParams.MinimalPowerdB=-40;        % Echos below this Limit are not initialised1

% ---- UserParameters ---

ChannelParams.DistanceFromRoadMiddle=-5; % negative: continental (right), positive: England (left)1

% ---- Graphics Parameters ---           

ChannelParams.GraphicalPlotArea=50;     % 1
ChannelParams.ViewVector = [-60,20];     % 3D visualization is not available in the free version1
ChannelParams.RoadWidth = 15;            %1

% --- Building Params ---

ChannelParams.BuildingRow1=1;            % logigal to switch Building Row right(heading 0 deg) on1
ChannelParams.BuildingRow2=1;            % logigal to switch Building Row left (heading 0 deg) on1
ChannelParams.BuildingRow1YPosition=-12; % m1
ChannelParams.BuildingRow2YPosition=12;  % m1

ChannelParams.HouseWidthMean=22;         % m1
ChannelParams.HouseWidthSigma=25;        % m1
ChannelParams.HouseWidthMin=10;          % m1
ChannelParams.HouseHeightMin=4;          % m1
ChannelParams.HouseHeightMax=50;         % m1
ChannelParams.HouseHeightMean=16;        % m1
ChannelParams.HouseHeightSigma=6.4;      % m1
ChannelParams.GapWidthMean=27;           % m1
ChannelParams.GapWidthSigma=25;          % m1
ChannelParams.GapWidthMin=10;            % m1
ChannelParams.BuildingGapLikelihood=0.18;% lin Value1

% --- Tree Params ---

ChannelParams.TreeHeight = 8;            % m
ChannelParams.TreeDiameter = 5;          % m
ChannelParams.TreeTrunkLength=2;         % m
ChannelParams.TreeTrunkDiameter=.2;      % m

ChannelParams.TreeAttenuation = 1.1;     % dB/m

ChannelParams.TreeRow1Use=1;             % logical switches tree row 1 on
ChannelParams.TreeRow2Use=1;             % logical switches tree row 2 on

ChannelParams.TreeRow1YPosition=-8;      % m
ChannelParams.TreeRow2YPosition=8;       % m

ChannelParams.TreeRow1YSigma=2;          % m
ChannelParams.TreeRow2YSigma=2;          % m

ChannelParams.TreeRow1MeanDistance=60;   % m
ChannelParams.TreeRow2MeanDistance=40;   % m

ChannelParams.TreeRow1DistanceSigma=20;  % m
ChannelParams.TreeRow2DistanceSigma=20;  % m

% --- Pole Params ---

ChannelParams.PoleHeight = 10;           % m
ChannelParams.PoleDiameter = .2;         % m

ChannelParams.PoleRow1Use=1;             % logical switches Pole row 1 on
ChannelParams.PoleRow2Use=0;             % logical switches Pole row 2 on

ChannelParams.PoleRow1YPosition=0;       % m
ChannelParams.PoleRow2YPosition=10;      % m

ChannelParams.PoleRow1YSigma=1;          % m
ChannelParams.PoleRow2YSigma=1;          % m

ChannelParams.PoleRow1MeanDistance=25;   % m
ChannelParams.PoleRow2MeanDistance=10;   % m

ChannelParams.PoleRow1DistanceSigma=10;  % m
ChannelParams.PoleRow2DistanceSigma=10;  % m

% ------------ Initial Settings -------------
% - Anything Below here must not be changed -
% -------------------------------------------

Co=2.99e8; % Speed of Light

MaximumPossibleSpeed=Co*Parameters.SampFreq/ChannelParams.CarrierFreq/2; % To fulfill the sampling Theorem
SamplingTime=1/Parameters.SampFreq;

% --- Initialising the channel object ---

pause(1)
disp('Initialising the channel ...')
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

h = waitbar(0,'Simulation running ...');

if ChannelParams.EnableCIRDisplay
    
    % --- init CIR figure ---

    hh = figure;
    subplot(211)
    xlabel('Delay in s')
    ylabel('Power in dB')
    axis([-2e-8,40e-8,0,50])
    set(get(hh,'Children'),'YTickLabel',[-40 -30 -20 -10 0 10]);
    hold on
    grid on
    plot (0,0,'r')
    plot (0,0)
    legend ('Direct paths','Echo paths')
    
    subplot(212)
    xlabel('Delay in s')
    ylabel('Phase in rad')
    axis([-2e-8,40e-8,-pi,pi])
    hold on
    grid on

end

for dhv=1:Parameters.NumberOfSteps

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
    waitbar(dhv/Parameters.NumberOfSteps,h)
    
    % --- binning LOS ---
    
    for sfg = 1:length(LOSDelays)
        
        dlybin = round(LOSDelays(sfg)/dlystp) + 1;
        pwrbin = round(20*log10(abs(LOS(sfg)))/pwrstp) + 1;
        
        if pwrbin<=length(pwrvec) & pwrbin>0 & dlybin<=length(dlyvec)
            PowerDelayProfile(pwrbin,dlybin) = PowerDelayProfile(pwrbin,dlybin) + 1;
        end
        
    end  
    
    % --- binning echoes ---
    
    for sfg = 1:length(DelayVec)
        
        dlybin = round(DelayVec(sfg)/dlystp) + 1;
        pwrbin = round(20*log10(abs(ComplexOutputVec(sfg)))/pwrstp) + 1;
        
        if pwrbin<=length(pwrvec) & pwrbin>0 & dlybin<=length(dlyvec)
            PowerDelayProfile(pwrbin,dlybin) = PowerDelayProfile(pwrbin,dlybin) + 1;
        end
        
    end  
    
    if ChannelParams.EnableCIRDisplay
        
        % --- display CIR ---

        figure(hh);
        subplot(211)
        cla
        Time = dhv/Parameters.SampFreq;
        title(['Channel Impulse Response (CIR), T = ',num2str(Time,'%5.2f'),' s, v = ',num2str(ActualSpeed*3.6,'%4.1f'),' km/h'])
        stem(LOSDelays,40 + 20*log10(abs(LOS)),'r');
        stem(DelayVec,40 + 20*log10(abs(ComplexOutputVec)));
        
        subplot(212)
        cla
        stem(LOSDelays,angle(LOS),'r');
        stem(DelayVec,angle(ComplexOutputVec));
        
    end

end%
close(h);

% --- calculate probability density function ---

PowerDelayProfile = PowerDelayProfile/sum(sum(PowerDelayProfile));

% --- display PowerDelayProfile ---

figure
surf(dlyvec,pwrvec,10*log10(PowerDelayProfile+eps),'LineStyle','none','FaceColor','interp','EdgeLighting','phong');
caxis([-70,-10]);
view(2)
xlabel('delay in s')
ylabel('power in dB')
title('Power delay profile - probability density function')
hc = colorbar;
set(hc,'YLim',[-70,-10])
clear newYTic YTic
YTic = get(hc,'YTickLabel');
% axes(hc);
ax = axis;
dta = (ax(3)-ax(4))/(size(YTic,1)-1);
for kk = 1:size(YTic,1)
    oldYTic(kk,:) = [' '];
    YTicString = YTic(kk,:);
    newYTic(kk,:) = ['10^{',num2str(str2num(YTicString{1})/10),'}'];
    % text(1.2,ax(3)-dta*(kk-1),['10^{',num2str(str2num(YTicString{1})/10),'}'],'interpreter','tex','horizontalalignment','left','verticalalignment','middle', FontSize=20);
end
set(hc,'YTickLabel',newYTic);

% ---------------------------------

disp(' ');
disp('Simulation finished');

