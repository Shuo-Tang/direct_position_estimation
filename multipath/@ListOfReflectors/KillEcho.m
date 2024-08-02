function [this,UnusedTime]=KillEcho(this,NumberOfEchosToBeKilled,CurrentRxPosition)

%
% Copyright © 2005
% Deutsches Zentrum für Luft- und Raumfahrt e.V.
% German Aerospace Centre
%
% Deletes N Echos of the list. Those are selected which lifetime is most
% expired.
% call: [this,UnusedTime]= KillEcho(this,NumberOfEchosToBeKilled,RxPosition)
%
% NumberOfEchosToBeKilled = Number of Echos to be killed
% CurrentRxPosition = Current receiver position
%
% UnusedTime = sum of remaining lifetime for Echos

EchoPosis=find(~this.EmptyCells);
TotalEchos = length(EchoPosis);

UnusedTime=0;

if NumberOfEchosToBeKilled >TotalEchos
    error('you cannot kill more Echos than are in the list')
end %if

if NumberOfEchosToBeKilled~=0
    for dhv=1:TotalEchos
        CurrentEcho(dhv)=EchoPosis(dhv);

        UsedLifetime(dhv)=norm(this.RxBirthCoordinate{CurrentEcho(dhv)}-CurrentRxPosition);
        RemainingLifeTime(dhv)=this.LifeTimeInM(CurrentEcho(dhv))-UsedLifetime(dhv);

    end % for

    [SortedRemainingLifetime,Index]=sort(RemainingLifeTime);

    EchosToBeKilled=Index(1:NumberOfEchosToBeKilled);
    IndexNumbersToBeKilled=EchoPosis(EchosToBeKilled);

    % Kill them

    this.EmptyCells(IndexNumbersToBeKilled)=1;
    UnusedTime=sum(max(0,SortedRemainingLifetime(1:NumberOfEchosToBeKilled)));
end %if




