function varargout=KnifeEdge(varargin)

%
% Copyright © 2005
% Deutsches Zentrum für Luft- und Raumfahrt e.V.
% German Aerospace Centre
%
% this function calculates knife-edge diffraction 
%
% call:
% 
%                       nu = KnifeEdge(d1,d2,h,h1,y1,theta,lambda);
%  [Amplitude,DetourLength]= KnifeEdge(d1,d2,h,h1,y1,theta,lambda);
%                Amplitude = KnifeEdge(nu);
%
% for details in notation see Example 16.4.1 in  
% 'Electromagnetic Waves and Antennas' by Sophocles J. Orfanidis,
% available at www.ece.rutgers.edu/~orfanidi/ewa
%

switch nargin
  case 1,
    nue=varargin{1};
  case 7
    d1=varargin{1}; d2=varargin{2}; h=varargin{3}; h1=varargin{4}; y1=varargin{5};
    theta=varargin{6}; lambda=varargin{7}; 
    F=d1.*d2/(d1+d2);  % focal length
    b1=y1.*cos(theta); % clearance distance
    alpha=b1./F;
    nue=alpha.*sqrt(2*F/lambda);
  otherwise
    error('Wrong number of input arguments.');
end;

if (nargout == 1) && (nargin == 7),
  varargout{1}=nue;
else
  if abs(nue)<1e9;
    %[C,S] = fcs(min(nue,1e9)); % Fresnel integral
    Z=fcs(min(nue,1e9)); % Fresnel integral
    C = real(Z);
    S = imag(Z);
    Amplitude=1/(1-i)*((C-i*S)+(1-i)/2); % diffraction coefficient (complex)
  else % outside valid range of fcs implementation
    Amplitude=1*(nue>0); % diffraction negligible 
  end;
  varargout{1}=Amplitude;
  if (nargout == 2) && (nargin == 7),
    DetourLength=sqrt((h-h1)^2+d1^2) - d1/cos(theta) + y1*sin(theta); % from geometry
    varargout{2}=DetourLength;
  elseif nargout > 2,
    error('Wrong number of output arguments.');
  end;
end;

