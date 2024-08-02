function NrOfEchos = GetNumberOfEchos(this)

%
% Copyright © 2005
% Deutsches Zentrum für Luft- und Raumfahrt e.V.
% German Aerospace Centre
%
% returns the number of Echos of the LIst
%
% call: 
%
% NrOfEchos = GetNumberOfEchos(Object)

NrOfEchos=this.MaximumNumberOfReflectors-sum(this.EmptyCells);

