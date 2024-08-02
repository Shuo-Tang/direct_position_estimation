function varargout=generate(this,ActualSpeed,ActualHeading,SatElevation,SatAzimut)

%
% Copyright © 2005
% Deutsches Zentrum für Luft- und Raumfahrt e.V.
% German Aerospace Center
%
% Generate routine of the Land Mobile Multipath Channel
%
% call:
%
% [this,LOS,LOSDelay,EchoValues,EchoDelays,EchoIdentifier,CurrentReceiverPosition,Time]=generate(this,ActualSpeed,ActualHeading,SatElevation,SatAzimut)
% 
% Please note: variable number of output Arguments!

Co=3e8;

% ---- calculating Way -----

WayForward=ActualSpeed/this.SampFreq;
this.CurrentRxPosition=this.CurrentRxPosition+[1 0 0]*WayForward;

% Conversion of Sat Azimut and Heading

RelativeBearing=ActualHeading-SatAzimut;

ElevRad=pi/2+pi*SatElevation/180;
AzimRad=pi+pi*RelativeBearing/180;

this.SatVector=[sin(ElevRad)*cos(AzimRad),sin(ElevRad)*sin(AzimRad),cos(ElevRad)];

% ---- calculating Direct path ----

ComplexValue(1)=1;
Delay(1)=0;
Identifier(1)=0;


% --- Checking the LifeTime of Echos

this.TheReflectorList=CheckLifetimeOfEchos(this.TheReflectorList,this.CurrentRxPosition);


% --- Calculating the current number of echos ----

[this.NumberGenerator,TargetNumberOfEchos] = generate(this.NumberGenerator,this.CurrentRxPosition,SatElevation);
  
% limit to maximum number of echoes
TargetNumberOfEchos = min(TargetNumberOfEchos,this.MaximumNumberOfReflectors);

ActualNumberOfEchos=GetNumberOfEchos(this.TheReflectorList);

EchoNumberDifference=ActualNumberOfEchos-TargetNumberOfEchos;


% ----- Determine if new Echos are necessary -----

if EchoNumberDifference < 0
    
    % ---- Determine Parameters of the new Echo ----

    EnvironmentParams.Azimuth           = RelativeBearing;
    EnvironmentParams.Elevation         = SatElevation;
    EnvironmentParams.BuildHeightMean   = this.HouseHeightMean;
    EnvironmentParams.BuildHeightSigma  = this.HouseHeightSigma;
    EnvironmentParams.RxPosition        = this.CurrentRxPosition;

    % create Echo Number Difference x New Echos

    for dhv=1:abs(EchoNumberDifference)

        [this.RandomEchoGenerator,ReflectorParams,LifeDistance,ReflectorSpeedRatio] = generate(this.RandomEchoGenerator,EnvironmentParams);

        % ---- Initialise the New Echo ----

        this.TheReflectorList=NewEcho(this.TheReflectorList,ReflectorParams,this.CurrentRxPosition,LifeDistance,ReflectorSpeedRatio);

    end%for
end%if

% ----- Determine if  Echos are to be killed -----

if EchoNumberDifference > 0
    this.TheReflectorList=KillEcho(this.TheReflectorList,EchoNumberDifference,this.CurrentRxPosition);
    %     disp('killing Echo')
end %if


% ---- calculation of the echos ----

[this.TheReflectorList,ComplexOutputVec,DelayVec,EchoNumberVec]=GenerateCIR(this.TheReflectorList,this.CurrentRxPosition,this.SatVector,this.Time);

% Direct path generation

LOS=1;
LOSDelay=[];

% ---- calculation of the trees ----

if this.TreeRow1Use
    this.TreeRow1=CreateNewKillOld(this.TreeRow1,this.CurrentRxPosition);
    LOS=LOS*generate(this.TreeRow1,this.CurrentRxPosition,this.SatVector,this.CarrierFreq);
end %if

if this.TreeRow2Use
    this.TreeRow2=CreateNewKillOld(this.TreeRow2,this.CurrentRxPosition);
    LOS=LOS*generate(this.TreeRow2,this.CurrentRxPosition,this.SatVector,this.CarrierFreq);
end %if
% ---- calculation of the Poles ----
if this.PoleRow1Use
    this.PoleRow1=CreateNewKillOld(this.PoleRow1,this.CurrentRxPosition);
    LOS=LOS*generate(this.PoleRow1,this.CurrentRxPosition,this.SatVector,this.CarrierFreq);
end% if

if this.PoleRow2Use
    this.PoleRow2=CreateNewKillOld(this.PoleRow2,this.CurrentRxPosition);
    LOS=LOS*generate(this.PoleRow2,this.CurrentRxPosition,this.SatVector,this.CarrierFreq);
end%if

% Please note: Since the user is inbetween Housefront 1 and the housefront
% 2 the LOS variables can be multiplied.
LOSDelay = 0;  % for the case that no buildings exist

if this.BuildingRow1
    this.HouseFront1=CreateNewKillOld(this.HouseFront1,this.CurrentRxPosition);
    [HouseFront1Amplitude,HouseFront1Delay]=generate(this.HouseFront1,this.CurrentRxPosition,this.SatVector,this.CarrierFreq);
    LOS=LOS*HouseFront1Amplitude;
    LOSDelay=HouseFront1Delay;
end%if

if this.BuildingRow2
    this.HouseFront2=CreateNewKillOld(this.HouseFront2,this.CurrentRxPosition);
    [Amplitudes,HouseFront2Delay]=generate(this.HouseFront2,this.CurrentRxPosition,this.SatVector,this.CarrierFreq);
    if ~this.BuildingRow1 | (isequal(HouseFront1Amplitude,1) && isequal(HouseFront1Delay,0))
      % Building row 2 only used if building row 1 is non-existing or gives
      % LOS with delay 0 and amplitude 1
        LOSDelay=HouseFront2Delay;
        LOS=LOS*Amplitudes;
    end
end%if


% calculating Doppler LOS

WayVec=[WayForward,0,0];

LOSRelativePhaseTurn=(this.SatVector*WayVec')*this.CarrierFreq/Co*2*pi;

this.LastLOSPhase=LOSRelativePhaseTurn+this.LastLOSPhase;

LOS=LOS*exp(j*this.LastLOSPhase);



% ---- calculating Time -----


this.Time=this.Time+1/this.SampFreq;

% remove those paths below minimal power

PathPowerdB=10*log10(abs(ComplexOutputVec).^2);
UsedPaths=PathPowerdB>=this.MinimalPowerdB;

LOSPathPowerdB=10*log10(abs(LOS).^2);
LOSUsedPaths=LOSPathPowerdB>=this.MinimalPowerdB;

varargout{1}=this;
varargout{2}=LOS(LOSUsedPaths);
varargout{3}=LOSDelay(LOSUsedPaths);
varargout{4}=ComplexOutputVec(UsedPaths);
varargout{5}=DelayVec(UsedPaths);
varargout{6}=EchoNumberVec(UsedPaths);
varargout{7}=this.CurrentRxPosition(1);
varargout{8}=this.Time;

