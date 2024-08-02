function this = CreateNewKillOld(this,CentralViewPoint)

%
% Copyright © 2005
% Deutsches Zentrum für Luft- und Raumfahrt e.V.
% German Aerospace Centre
%
% update contents of scenery according to new position
%

CurrentPosition=CentralViewPoint(1);

% remove old values
ind=find(this.XValues<CurrentPosition-this.ScenerySize);
if ~isempty(ind),
  this.XValues(ind)=[];
  this.HeightValues(ind)=[];
  this.IDValuesX(ind)=[];
  this.IDValuesZ(ind)=[];
end;

% fill up scenery width with new positions
if isempty(this.XValues),
  this.XValues=CurrentPosition-this.ScenerySize;
  this.HeightValues=0;
  this=CreateIDs(this);
end;

while this.XValues(end) < CurrentPosition + this.ScenerySize, % generate new value

  if rand < this.GapLikelihood, % gap
    
    width=-1;
    while width < this.GapWidthMin,
      width=this.GapWidthMean+randn*this.GapWidthSigma;
    end;
    this.XValues(end+1)=this.XValues(end)+width;
    this.HeightValues(end+1)=0;
    this=CreateIDs(this);

  else % house

    width=-1;
    while width < this.HouseWidthMin,
      width=this.HouseWidthMean+randn*this.HouseWidthSigma;
    end;
    this.XValues(end+1)=this.XValues(end)+width;

    height=-1;
    while height < this.HouseHeightMin,
      height=min(this.HouseHeightMax,this.HouseHeightMean+randn*this.HouseHeightSigma);
    end;
    this.HeightValues(end+1)=height;
    this=CreateIDs(this);
  end;
end;

% helper function

function this=CreateIDs(this)

this.IDValuesX(end+1)=this.IDCounter;
this.IDValuesZ(end+1)=this.IDCounter+1;
this.IDCounter=this.IDCounter+2;

  

  
