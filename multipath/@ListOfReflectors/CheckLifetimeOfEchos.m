function [this]=CheckLifeTimeOfEchos(this,CurrentRxPosition)

%
% Copyright © 2005
% Deutsches Zentrum für Luft- und Raumfahrt e.V.
% German Aerospace Centre
%
% Deletes those Echos which LifeTime has expired
%
% call:
%
% [this]=CheckLifeTimeOfEchos(this,CurrentRxPosition)


EchoPosis=find(~this.EmptyCells);
TotalEchos = length(EchoPosis);


if ~all(this.EmptyCells)

    for dhv=1:TotalEchos
        CurrentEcho(dhv)=EchoPosis(dhv);

        UsedLifetime(dhv)=norm(this.RxBirthCoordinate{CurrentEcho(dhv)}-CurrentRxPosition);
        RemainingLifeTime(dhv)=this.LifeTimeInM(CurrentEcho(dhv))-UsedLifetime(dhv);

    end % for

    EchosToBeKilled=RemainingLifeTime<=0;
    IndexNumbersToBeKilled=EchoPosis(EchosToBeKilled);

    % Kill them

    this.EmptyCells(IndexNumbersToBeKilled)=1;

end %if