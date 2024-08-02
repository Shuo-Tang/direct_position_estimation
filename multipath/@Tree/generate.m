% this function generates a complex amplitude
%
% call:
%
% Amplitude=generate(TreeObject,RxVec,SatVec,CarrierFreq)
%


function Amplitude=generate(this,RxVec,SatVec,CarrierFreq)

HeightAtTree=[];
ComplexValueStatistical=1;

SatDistance=2e7; % distance to satellite
Lambda=3e8/CarrierFreq;

% consider projection on z=RxZ plane
SatVec=-SatVec/norm(SatVec); % unity vector from Rx to satellite
AziSat=atan2(SatVec(2),SatVec(1));
TreeRxVec=this.Position(1:2)-RxVec(1:2); % from Rx to tree on plane
AziTree=atan2(TreeRxVec(2),TreeRxVec(1));
alpha=AziSat-AziTree;

if abs(alpha) >= pi/2, % tree on other side as satellite?
  ADiffraction = 1;
  ADamping = 1;
else
  TreeDistance=norm(TreeRxVec);
  d=TreeDistance*sin(alpha); % closest distance between LOS and tree

  % check height
  SatVecProjection=TreeDistance*cos(alpha);
  Elevation=asin(SatVec(3));
  HeightAtTree=SatVecProjection*tan(Elevation)+RxVec(3);

  % diffraction at tree trunk?
  if (HeightAtTree < this.Position(3)) || (HeightAtTree > this.Position(3)+this.TreeTrunkLength)
    ADiffraction = 1; % LOS below or above tree trunk
  else
    Radius=0.5*this.TreeTrunkDiameter;
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

  % effects at tree top
  if (abs(d) > 0.5*this.Diameter),
    ADamping=1; % LOS does not touch
  else
    % height of chord within cylinder
    ChordProjection=sqrt((0.5*this.Diameter)^2-d^2);
    HeightOfChord=ChordProjection*tan(Elevation);
    HTop=this.Position(3)+this.ObstacleHeight;
    HMid=this.Position(3)+this.TreeTrunkLength;
    ChordTop=HeightAtTree+HeightOfChord;
    ChordBot=HeightAtTree-HeightOfChord;

    % determine part within tree top
    if (HTop < ChordBot) || (HMid > ChordTop),
      ADamping = 1; % no overlap
    else
      % damping at tree top
      if Elevation == 0,
        DampingdB=ChordProjection*this.MeanDampingFactor;
      else
        HeightInTop=min(HTop,ChordTop)-max(HMid,ChordBot);
        DampingdB=HeightInTop/sin(Elevation)*this.MeanDampingFactor;
      end;
      ADamping=10^(-DampingdB/10);

      % stochastic process at tree top
      ComplexValueStatistical = this.ProcessScalingFactor*(this.MyRiceanConst + (sqrt(2)*sum(exp(j*(this.RandomFrequencies*RxVec(1)*2*pi+this.RandomPhases))) / sqrt(length(this.RandomFrequencies))));
    end;
  end;
end;

% combine effect from trunk and top
Amplitude=ADiffraction*ADamping*ComplexValueStatistical;



