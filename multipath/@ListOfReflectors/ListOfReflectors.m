function obout = ListOfReflectors(Params)

%
% Copyright © 2005
% Deutsches Zentrum für Luft- und Raumfahrt e.V.
% German Aerospace Centre
%
% Constructor of List Of Reflectors
%
% call: Mylist=ListOfReflectors(Params);
%
% Params.MaximumNumberOfReflectors  = Gives the Maximum Number of Reflectors
% Params.EnableDisplay              = 1 Enables DIsplay for all Reflectors
% Params.CarrierFreq                = CarrierFrequency in Hz


% Please note: the only way an Echo is deleted: set EmptyCells=0. Then the
% Echo is not used anymore. 

new=Params;

% The List of Reflectors
new.Reflectors{Params.MaximumNumberOfReflectors}=[];

% List of Bith Coordinates of the Reflectors
new.RxBirthCoordinate{Params.MaximumNumberOfReflectors}=[];

% List of Livetime of the reflector in meters
new.LifeTimeInM=zeros(1,new.MaximumNumberOfReflectors);

% List of Empty Cells
% Value is one if equivalent Refletor is empty

new.EmptyCells=ones(1,new.MaximumNumberOfReflectors);

% List of EchoNumbers to identyfy the Echo later on

new.EchoNumbers=zeros(1,new.MaximumNumberOfReflectors);

% List of Reflector Speed ratios

new.ReflectorSpeedRatio=zeros(1,new.MaximumNumberOfReflectors);

% to be able to label an echo the last given number must be stored
new.LastEchoNumerUsed=0;

obout=class(new,'ListOfReflectors');