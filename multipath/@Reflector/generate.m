function [this,ComplexValue,Delay,DeterministicAbsolutePhase]=generate(this,RxPosition,SatVector,time,ReflectorMovement);

%
% Copyright © 2005
% Deutsches Zentrum für Luft- und Raumfahrt e.V.
% German Aerospace Centre
%
% method generate of class reflector
%
% generates the new complex sample in dependence of the new position
%
% function call: [MyEcho,CompexValue,Delay]=generate(MyEcho,RxPosition,SatVector,time)
%
% Argument: Receiver position [x y z] real vector in m
%           SatVector is a [x y z] vector pointing from the receiver to the
%           satellite
%           Speed in m/s
%           Time in sec
%
% Return: ComplexValue for the Delay in sec, DeterministicAbsolutePhase in rad 
%

% calculating the movement way of the reflector. 

ReflectorMovementVec=[ReflectorMovement,0,0];

 % General Calculations

 Co=2.99e8;
 
 % ----------------
 %      Delay
 % ----------------
 
 % calculating reflector to receiver vector
 RfxRxVec = RxPosition-(this.Position+ReflectorMovementVec);

 % calculating path delay difference of parallel wavefront
 % (projection of RfxRxVec in SatVector)
 Lambda = RfxRxVec*(-SatVector)';

 % path delay in m
 Delay = Lambda + norm(RfxRxVec);

 Delay=Delay/Co;

 
 % --------------------------------
 %      Complex Value
 % --------------------------------
 
 % --- statistical process ---


 ComplexValueStatistical = this.ProcessScalingFactor*(this.MyRiceanConst + (sqrt(2)*sum(exp(j*(this.RandomFrequencies*time*2*pi+this.RandomPhases))) / sqrt(length(this.RandomFrequencies))));
 
 
 % --- deterministic part ---
 
 
 BaselineRx_Refl=((this.Position+ReflectorMovementVec)-RxPosition);
 Distance=norm(BaselineRx_Refl);
 DopplerPhase=norm(BaselineRx_Refl)/this.Wavelength*2*pi; 

 % --- Merge Deterministic, statistic --------
 
ComplexValue=exp(j*(-DopplerPhase+this.Phase))*ComplexValueStatistical;

% ----- Output of absolute Phase ----

DeterministicAbsolutePhase=-DopplerPhase+this.Phase;

this.ActualPowerdB=10*log10(abs(ComplexValue)^2);

