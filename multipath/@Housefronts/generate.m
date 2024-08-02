% this function generates a vector of complex amplitudes and a vector of Delays.
%
% call:
%
% [Amplitudes,PathDelays]=generate(HousfrontObject,RxVec,SatVec,CarrierFreq)
% [Amplitudes,PathDelays,Identifiers]=generate(HousfrontObject,RxVec,SatVec,CarrierFreq)
%

function varargout=generate(this,RxVec,SatVec,CarrierFreq)

% default: only LOS directly
Amplitudes=1;
PathDelays=0;
Identifiers=0;

% satellite on same side as houses?
if ~isequal(SatVec(2),0) && (sign(-SatVec(2)) == sign(this.YPosition-RxVec(2))), 

  SatDistance=2e7; % distance to satellite
  Lambda=3e8/CarrierFreq;

  % LOS intersection with housefront
  SatVec=-SatVec/norm(SatVec); % unity vector from Rx to satellite
  IntersectDistance=(this.YPosition-RxVec(2))/SatVec(2); % from Rx to intersection point
  XAtHouse=IntersectDistance*SatVec(1)+RxVec(1);
  ZDistance=IntersectDistance*SatVec(3);
  ZAtHouse=ZDistance+RxVec(3);
  Elevation=asin(SatVec(3)); % elevation
  HouseDistance=IntersectDistance*cos(Elevation); % on ground plane
  YDistance=abs(this.YPosition-RxVec(2));

  AziStart=atan2(YDistance,this.XValues(end)-RxVec(1));
  AziEnd=atan2(YDistance,this.XValues(1)-RxVec(1));
  Azimuth=atan2(abs(SatVec(2)),SatVec(1));
  if (Azimuth > AziStart) && (Azimuth < AziEnd),  % satellite behind housefront

    % heights and walls of house around LOS
    ind=find(this.XValues >= XAtHouse); ind=ind(1);
    Height=this.HeightValues(ind);
    XLeft=-Inf; XRight=Inf;

    if ZAtHouse > Height,  % LOS does not touch houses

      % top of house below
      nue(1)=KnifeEdge(HouseDistance,SatDistance,Height,RxVec(3),ZAtHouse-Height,Elevation,Lambda);
      ID(1)=this.IDValuesZ(ind);

      % next wall to the left
      indl=find(this.HeightValues(1:ind-1) > ZAtHouse);
      d1=sqrt(YDistance^2+ZDistance^2);
      if ~isempty(indl),
        XLeft=this.XValues(indl(end));
        nue(2)=KnifeEdge(d1,SatDistance,XLeft,RxVec(1),XAtHouse-XLeft,atan((XAtHouse-RxVec(1))/d1),Lambda);
        ID(2)=this.IDValuesX(indl(end));
      else
        nue(2)=Inf;
        ID(2)=-Inf;
      end;

      % next wall to the right
      indr=find(this.HeightValues(ind+1:end) > ZAtHouse);
      if ~isempty(indr),
        XRight=this.XValues(ind+indr(1)-1);
        nue(3)=KnifeEdge(d1,SatDistance,-XRight,-RxVec(1),XRight-XAtHouse,atan((RxVec(1)-XAtHouse)/d1),Lambda);
        ID(3)=this.IDValuesX(ind+indr(1)-1);
      else
        nue(3)=Inf;
        ID(3)=Inf;
      end;

      % choose relevant edge for diffraction
      [dummy,closest]=min(nue);
      if (Height==0) && (closest==1),
        Amplitudes(1)=1;  % inside a gap
        Identifiers(1)=0; % LOS path directly
      else
        Identifiers(1)=ID(closest);
        Amplitudes(1)=KnifeEdge(nue(closest));
      end;
      PathDelays(1)=0; % no delay since LOS is there

    else % LOS traverses through house

      % top of house
      [Amplitudes(1),DetourLength]=KnifeEdge(HouseDistance,SatDistance,Height,RxVec(3),ZAtHouse-Height,Elevation,Lambda);
      PathDelays(1)=DetourLength/3e8;
      Identifiers(1)=this.IDValuesZ(ind);

      % next wall to the left
      d1=sqrt(YDistance^2+ZDistance^2);
      indl=find(this.HeightValues(1:ind-1) < ZAtHouse);
      if ~isempty(indl),
        XLeft=this.XValues(indl(end));
        [Amplitudes(2),DetourLength]=KnifeEdge(d1,SatDistance,-XLeft,-RxVec(1),XLeft-XAtHouse,atan((RxVec(1)-XAtHouse)/d1),Lambda);
        PathDelays(2)=DetourLength/3e8;
        Identifiers(2)=this.IDValuesX(indl(end));
      end;

      % next wall to the right
      indr=find(this.HeightValues(ind+1:end) < ZAtHouse);
      if ~isempty(indr),
        XRight=this.XValues(ind+indr(1)-1);
        [Amplitudes(3),DetourLength]=KnifeEdge(d1,SatDistance,XRight,RxVec(1),XAtHouse-XRight,atan((XAtHouse-RxVec(1))/d1),Lambda);
        PathDelays(3)=DetourLength/3e8;
        Identifiers(3)=this.IDValuesX(ind+indr(1)-1);
      end;

    end;
  end;
end;

varargout{1}=Amplitudes;
varargout{2}=PathDelays;
varargout{3}=Identifiers;


