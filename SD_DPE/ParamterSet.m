clear
close all
clc
%% Set SVs' position and velocity
numSV = 7;
GNSS_config.no_sat = numSV+1;
time = 0;
[SatPosition,SatVelocity] = Satellite_positions_and_velocities(time,...
    GNSS_config);
%% Set Clock bias for satellites
deltat =  randn(1,numSV)*1e-10;
%% Set user's position and velocity (for instance: 42?20'13.4"N
% 71?05'25.9"W)
true_phi_b_dms = [42 20 13.4];        % latitude (degree, minutes, seconds)
true_lambda_b_dms = -[71 05 25.9];    % longitude (degree, minutes, seconds)
true_h_b = 10;                        % height (meters)
% some transformations
true_phi_b_deg = true_phi_b_dms(1) + true_phi_b_dms(2)/60 + true_phi_b_dms(3)/3600;% latitude (decimal degrees)
true_phi_b_rad = pi/180*true_phi_b_deg;% latitude (radians)

true_lambda_b_deg = true_lambda_b_dms(1) + true_lambda_b_dms(2)/60 + true_lambda_b_dms(3)/3600;% longitude (decimal degrees)
true_lambda_b_rad = pi/180*true_lambda_b_deg;% longitude (radians)

UserPosition = pv_NED_to_ECEF(true_phi_b_rad,true_lambda_b_rad,true_h_b,[0;0;0])';% Ellipsoidal-to-Cartesian
UserVelocity = [0 0 0]; % suppose that target is static
%% Set Clock bias and drift for the user
deltaT = randn*1e-6;
deltaTdot = abs(randn*1e-10); % 0.1 nano drift in one second
%% Compute Unit Direction
UnitDirection = (SatPosition - UserPosition)./vecnorm(SatPosition - UserPosition,2,2); % unit vector of direction
%% save parameters
param.SatPosition = SatPosition(2:numSV+1,:);
param.SatVelocity = SatVelocity(2:numSV+1,:);
param.deltat = deltat;
param.UserPosition = UserPosition;
param.UserVelocity = UserVelocity;
param.deltaT = deltaT;
param.deltaTdot = deltaTdot;
param.UnitDirection = UnitDirection(2:numSV+1,:);
save('param.mat','param')