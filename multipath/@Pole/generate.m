function ADiffraction=generate(this,RxVec,SatVec,CarrierFreq)

%
% Copyright © 2005
% Deutsches Zentrum für Luft- und Raumfahrt e.V.
% German Aerospace Centre
%
% this function generates a complex amplitude
%
% call:
%
% Amplitude=generate(PoleObject,RxVec,SatVec,CarrierFreq)
%

SatDistance=2e7; % distance to satellite
Lambda=3e8/CarrierFreq;

% consider projection on z=RxZ plane
SatVec=-SatVec/norm(SatVec); % unity vector from Rx to satellite
AziSat=atan2(SatVec(2),SatVec(1));
PoleRxVec=this.Position(1:2)-RxVec(1:2); % from Rx to pole on plane
AziPole=atan2(PoleRxVec(2),PoleRxVec(1));
alpha=AziSat-AziPole;

PoleDistance=norm(PoleRxVec);
d=PoleDistance*sin(alpha); % closest distance between LOS and pole

if abs(alpha) >= pi/2,
  ADiffraction = 1;  % too far away for diffraction
else
  % check height
  SatVecProjection=PoleDistance*cos(alpha);
  Elevation=asin(SatVec(3));
  HeightAtPole=SatVecProjection*tan(Elevation)+RxVec(3);

  if (HeightAtPole < this.Position(3)) || (HeightAtPole > this.Position(3)+this.ObstacleHeight)
    ADiffraction = 1; % LOS below or above pole
  else
    Radius=0.5*this.Diameter;
    d1=SatVecProjection/cos(Elevation);
    F=d1*SatDistance/(d1+SatDistance); % focal length
    Scaling=sqrt(2/(Lambda*F));
    
    % ensure smooth transitions 
    if ~isempty(this.FresnelRange), % use a look-up table
      IndR=round((Radius*Scaling-this.FresnelRange.MinR)/this.FresnelRange.RSpacing+1); 
      try
        NueMax=this.FresnelRange.Table(IndR);
      catch
        NueMax=FresnelRange(Radius*Scaling);
      end;
    else % search for a root
      NueMax=FresnelRange(Radius*Scaling);
    end;
    
    if abs(d*Scaling) < NueMax, % inside 
      NueLeft=(d-Radius)*Scaling;
      NueRight=(-Radius-d)*Scaling;
      ADiffraction=KnifeEdge(NueLeft)+KnifeEdge(NueRight);
    else
      ADiffraction=1;
    end;

  end;
end;





