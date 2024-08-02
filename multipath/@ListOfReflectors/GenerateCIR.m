function varargout=GenerateCIR(this,RxPosition,SatVector,time)

%
% Copyright © 2005
% Deutsches Zentrum für Luft- und Raumfahrt e.V.
% German Aerospace Centre
%
% this function generates a vector of Complex values and a vector of Delays.
%
% call:
%
% [MyList,ComplexVector,PathDelays]=GenerateCIR(MyList,RxPosition,SatVector,time)
% [MyList,ComplexVector,PathDelays,EchoNumbers]=GenerateCIR(MyList,RxPosition,SatVector,time)
% 

DelayVec=[];
ComplexOutputVec=[];

ActiveEchoIndexes=find(~this.EmptyCells);

for dhv=1:length(ActiveEchoIndexes)
    % calculates the Movement of the reflector due to is given speed ratio.
    % 
    
    
    ReflectorMovement=(RxPosition(1)-this.RxBirthCoordinate{ActiveEchoIndexes(dhv)}*[1 0 0]')*this.ReflectorSpeedRatio(ActiveEchoIndexes(dhv));
    [this.Reflectors{ActiveEchoIndexes(dhv)},ComplexOutputVec(dhv),DelayVec(dhv),DeterministicAbsolutePhase]=generate(this.Reflectors{ActiveEchoIndexes(dhv)},RxPosition,SatVector,time,ReflectorMovement);
end %for

varargout{1}=this;
varargout{2}=ComplexOutputVec;
varargout{3}=DelayVec;
varargout{4}=this.EchoNumbers(ActiveEchoIndexes);
this.EchoNumbers(ActiveEchoIndexes);
