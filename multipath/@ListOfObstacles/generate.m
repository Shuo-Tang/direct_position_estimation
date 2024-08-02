function out=generate(this,RxVec,SatVec,CarrierFreq)

%
% Copyright © 2005
% Deutsches Zentrum für Luft- und Raumfahrt e.V.
% German Aerospace Centre
%
% Generates the Components of the List of obstacles
%

out=1;

for dhv=1:length(this.TheList)
        
    out=out*generate(this.TheList{dhv},RxVec,SatVec,CarrierFreq);
    
end %for

