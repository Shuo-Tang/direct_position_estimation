function argout=FresnelRange(varargin)

%
% Copyright © 2005
% Deutsches Zentrum für Luft- und Raumfahrt e.V.
% German Aerospace Centre
%
% Computes diffraction range around poles and tree trunks.
% For smooth transitions, NueMax is chosen as root of Fresnel function. 
%
% Variant 1:    NueMax = FresnelRange(Radius);
%
% Variant 2:    NueMaxTable = FresnelRange(MinR,MaxR,Spacing);
%
%   with  MinR      : smallest Radius allowed
%         MaxR      : largest Radius allowed
%         RSpacing  : spacing

% rough location of desired roots is given by following table
NueLimits=[ ...,        
  0.1000   10.0000;  0.1100    9.1000;  0.1200    8.3000;  0.1300    7.7000;  0.1400    7.1000; ...
  0.1500    6.7000;  0.1600    6.2000;  0.1700    5.8000;  0.1800    5.5000;  0.1900    5.1000; ...
  0.2000   10.1000;  0.2100    9.7000;  0.2200    9.2000;  0.2300    8.8000;  0.2400    8.4000; ...
  0.2500    8.0000;  0.2600    7.8000;  0.2700    7.6000;  0.2800    7.3000;  0.2900    7.0000; ...
  0.3000   10.0000;  0.3100    9.7000;  0.3200    9.5000;  0.3300    9.2000;  0.3400    8.9000; ...
  0.3500    8.6000;  0.3600    8.4000;  0.3700    8.2000;  0.3800    8.1000;  0.3900    7.8000; ...
  0.4000   10.0000;  0.4100    9.9000;  0.4200    9.6000;  0.4300    9.4000;  0.4400    9.2000; ...
  0.4500    9.0000;  0.4600    8.8000;  0.4700    8.7000;  0.4800    8.4000;  0.4900    8.3000; ...
  0.5000   10.0000;  0.5100    9.8000;  0.5200    9.7000;  0.5300    9.5000;  0.5400    9.3000; ...
  0.5500    9.2000;  0.5600    9.0000;  0.5700    8.9000;  0.5800    8.8000;  0.5900    8.5000; ...
  0.6000   10.0000;  0.6100    9.9000;  0.6200    9.7000;  0.6300    9.6000;  0.6400    9.4000; ...
  0.6500    9.3000;  0.6600    9.2000;  0.6700    9.1000;  0.6800    8.9000;  0.6900    8.7000; ...
  0.7000   10.0000;  0.7100    9.9000;  0.7200    9.8000;  0.7300    9.6000;  0.7400    9.5000; ...
  0.7500    9.4000;  0.7600    9.3000;  0.7700    9.2000;  0.7800    9.1000;  0.7900    9.0000; ...
  0.8000   10.0000;  0.8100    9.9000;  0.8200    9.8000;  0.8300    9.7000;  0.8400    9.6000; ...
  0.8500    9.5000;  0.8600    9.4000;  0.8700    9.3000;  0.8800    9.2000;  0.8900    9.1000; ...
  0.9000   10.0000;  0.9100    9.9000;  0.9200    9.8000;  0.9300    9.7000;  0.9400    9.6000; ...
  0.9500    9.5000;  0.9600    9.4000;  0.9700    9.3000;  0.9800    9.2000;  0.9900    9.1000; ...
  1.0000   10.0000;
  ];

switch nargin
  case 1 % search root for given Radius

    [dummy,pos]=min(abs(NueLimits(:,1)-varargin{1}));
    argout=fzero(@RootsTwoSidedFresnel,NueLimits(pos,2),optimset('fzero'),varargin{1});
    
  case 3, % create look-up table

    Radius=varargin{1}:varargin{3}:varargin{2}; % grid for roots
    IndFine=find((Radius<=1) & (Radius >=0.1)); % range with specified root locations
    N=length(Radius);
    Limit=zeros(1,N);

    % interpolate and find closest root
    StartPos=10*ones(1,N); % default value
    StartPos(IndFine) = interp1(NueLimits(:,1),NueLimits(:,2),Radius(IndFine));    
    for i=1:N,
      Limit(i)=fzero(@RootsTwoSidedFresnel,StartPos(i),optimset('fzero'),Radius(i));
    end;
    argout.MinR=varargin{1};
    argout.MaxR=varargin{2};
    argout.RSpacing=varargin{3};
    argout.Table=Limit;

  otherwise
    error('Wrong number of input arguments.');

end;



% helper function
function y=RootsTwoSidedFresnel(Nue,Radius)

y=abs(KnifeEdge(Nue-Radius)+KnifeEdge(-Nue-Radius))-1;






