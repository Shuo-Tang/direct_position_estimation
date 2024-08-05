clear
% close all
% clc
rng(100)
%% Set SVs' position and velocity
numSV = 100;
GNSS_config.no_sat = numSV+1;
time = 0;
[SatPosition,SatVelocity] = Satellite_positions_and_velocities(time,...
    GNSS_config);
%% Set Clock bias for satellites
deltat =  randn(1,numSV+1)*1e-10;
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
%% Elevation Mask Check
true_p_eb_ecef = pv_NED_to_ECEF(true_phi_b_rad,true_lambda_b_rad,true_h_b,[0;0;0]);                     % Ellipsoidal-to-Cartesian
% Mask angle (deg)
GNSS_config.mask_angle = 10;
cos_lat = cos(true_phi_b_rad);
sin_lat = sin(true_phi_b_rad);
cos_long = cos(true_lambda_b_rad);
sin_long = sin(true_lambda_b_rad);
C_e_n = [-sin_lat * cos_long, -sin_lat * sin_long,  cos_lat;...
                   -sin_long,            cos_long,        0;...
         -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat];
elevation = zeros(numSV, 1);
for i = 1:numSV+1
    % Determine ECEF line-of-sight vector using (8.41)
    delta_r = SatPosition(i, :)' - true_p_eb_ecef;
    approx_range = sqrt(delta_r' * delta_r);
    u_as_e = delta_r / approx_range;            % line-of-sight unit vector
    
    % Convert line-of-sight vector to NED using (8.39) and determine
    % elevation using (8.57)
    elevation(i, 1) = -asin(C_e_n(3,:) * u_as_e);
end
ind1 = find(elevation >= degtorad(GNSS_config.mask_angle));
ind2 = find(SatPosition(:,1).*SatPosition(:,2).*SatPosition(:,3)>0);
ind = intersect(ind1, ind2);
ind = ind(1:2:end);
%% save parameters
param.SatPosition = SatPosition(ind,:);
param.SatVelocity = SatVelocity(ind,:);
param.deltat = deltat(:, ind);
param.UserPosition = UserPosition;
param.UserVelocity = UserVelocity;
param.deltaT = deltaT;
param.deltaTdot = deltaTdot;
param.UnitDirection = UnitDirection(ind,:);
save('param.mat','param')