function this=NewEcho(this,Params,RxPosition,LifeTimeInM,ReflectorSpeedRatio)

%
% Copyright © 2005
% Deutsches Zentrum für Luft- und Raumfahrt e.V.
% German Aerospace Centre
%
% creates a new reflector for given Parameters and includes it into a list.
%
% call:
%
% MyList=NewEcho(MyList,Params,RxPosition,LifeTimeInM,ReflectorSpeedRatio)
% 
% Params = Parameters of a Reflector object
% Rx Position = Rx Position at bith of Reflector
% LifeTimeInM = Lifetime (way) of the reflector in m
% ReflectorSpeedRatio  = the ratio between the receiver speed and the speed of the Reflector. Usually 0 or 1


if(all(this.EmptyCells==0))
    error('ListOfReflectors is full pleace watch the number of already existing reflectors next time')
end % if

PossiblePositions=find(this.EmptyCells==1);
ThePosition=PossiblePositions(1);

this.EmptyCells(ThePosition)=0;

Params.EnableDisplay=this.EnableDisplay;
Params.CarrierFreq=this.CarrierFreq;


this.Reflectors{ThePosition}=Reflector(Params);
this.RxBirthCoordinate{ThePosition}=RxPosition;
this.LifeTimeInM(ThePosition)=LifeTimeInM;
this.ReflectorSpeedRatio(ThePosition)=ReflectorSpeedRatio;

% set Echo Number

this.EchoNumbers(ThePosition)=this.LastEchoNumerUsed+1;
this.LastEchoNumerUsed=this.LastEchoNumerUsed+1;