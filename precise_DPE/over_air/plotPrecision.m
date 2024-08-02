% Javier Arribas 2012
% Accuracy and Precision analysis
%close all;
figure;
clear E;
clear N;
clear U;
% 1. Convert the GNSS Cartesian coordinates ECEF (X,Y,Z) to local topocentric
% coordinates (East,North,Up)

mean_Latitude=mean(navSolutions.latitude);
mean_Longitude=mean(navSolutions.longitude);
utmZone = findUtmZone(mean_Latitude,mean_Longitude);

% Compute moving averages because the recorded PVT is RAW data
averaging_depth=1000;

MX=tsmovavg(navSolutions.X, 's', averaging_depth);
MY=tsmovavg(navSolutions.Y, 's', averaging_depth);
MZ=tsmovavg(navSolutions.Z, 's', averaging_depth);

numPoints=length(navSolutions.X);
start_offset=1000;
step=1;
end_offset=start_offset+30000;
%end_offset=numPoints;
k=0;
for n=start_offset:step:end_offset
    k=k+1;
    %[E(k),N(k),U(k)]=cart2utm(navSolutions.X(n), navSolutions.Y(n), navSolutions.Z(n), utmZone);
    [E(k),N(k),U(k)]=cart2utm(MX(n), MY(n), MZ(n), utmZone);
end

v_2d=[E;N].'; %2D East Nort position vectors
v_3d=[E;N;U].'; %2D East Nort position vectors

% 2D analysis
% Simulated X,Y measurements
%v1=randn(1000,2);

% 2D Mean and Variance
mean_2d=[mean(v_2d(:,1)) ; mean(v_2d(:,2))];
sigma_2d=[sqrt(var(v_2d(:,1))) ; sqrt(var(v_2d(:,2)))];

sigma_ratio_2d=sigma_2d(2)/sigma_2d(1);

% if sigma_ratio=1 -> Prob in circle with r=DRMS -> 65% 
DRMS=sqrt(sigma_2d(1)^2+sigma_2d(2)^2);
% if sigma_ratio=1 -> Prob in circle with r=2DRMS -> 95% 
TWO_DRMS=2*DRMS;
% if sigma_ratio>0.3 -> Prob in circle with r=CEP -> 50%
CEP=0.62*sigma_2d(1)+0.56*sigma_2d(2);

scatter(v_2d(:,1)-mean_2d(1),v_2d(:,2)-mean_2d(2),'x');
hold on;

plot(0,0,'k*');

[x,y,z] = cylinder(DRMS,200);
plot(x(1,:),y(1,:),'r');
str = {'DRMS'};
text(cosd(45)*DRMS,cosd(45)*DRMS,str);

[x,y,z] = cylinder(TWO_DRMS,200);
plot(x(1,:),y(1,:),'g');
str = {'2DRMS'};
text(cosd(45)*TWO_DRMS,cosd(45)*TWO_DRMS,str);

[x,y,z] = cylinder(CEP,200);
plot(x(1,:),y(1,:),'r--');
str = {'CEP'};
text(cosd(45)*CEP,cosd(45)*CEP,str);

axis equal;

% 3D analysis
% Simulated X,Y,Z measurements
% v_3d=randn(1000,3);

% Mean and Variance
mean_3d=[mean(v_3d(:,1)) ; mean(v_3d(:,2)) ; mean(v_3d(:,3))];
sigma_3d=[sqrt(var(v_3d(:,1))) ; sqrt(var(v_3d(:,2))) ; sqrt(var(v_3d(:,3)))];


% if sigma_ratio=1 -> Prob in circle with r=DRMS -> 65% 
DRMS_3D=sqrt(sigma_3d(1)^2+sigma_3d(2)^2++sigma_3d(3)^2);
% if sigma_ratio=1 -> Prob in circle with r=2DRMS -> 95% 
TWO_DRMS_3D=2*DRMS_3D;

figure;
scatter3(v_3d(:,1)-mean_3d(1),v_3d(:,2)-mean_3d(2), v_3d(:,3)-mean_3d(3), 'x');
hold on;

%scatter3(v_3d(:,1)-mean_3d(1),v_3d(:,2)-mean_3d(2), -5*ones(numPoints,1), 'x','MarkerEdgeColor',[0.7 0.7 0.7]);
%scatter3(v_3d(:,1)-mean_3d(1),5*ones(numPoints,1), v_3d(:,3)-mean_3d(3), 'x','MarkerEdgeColor',[0.7 0.7 0.7]);

[x,y,z] = sphere();
hSurface=surf(DRMS_3D*x,DRMS_3D*y,DRMS_3D*z);  % sphere centered at origin
set(hSurface,'FaceColor',[0 1 0],'EdgeColor',[0 0.6 0],'EdgeAlpha',0.5,'FaceAlpha',0.4);

axis equal;

% Transform ECEF DOPS (the natural DPOP measurement obtained with from PVT) to local
% topocentric ENUP DOPS

ref_long=mean_Longitude;
ref_lat=mean_Latitude;

navSolutionsA_ECEF=mean(navSolutions.Q_DOP,3);

A_ECEF=A_ECEF(1:end-1,1:end-1); % exclude the time covariance

F=[-sind(ref_long) -sind(ref_lat)*cosd(ref_long) cosd(ref_lat)*cosd(ref_long); ...
  cosd(ref_long) -sind(ref_lat)*sind(ref_long) cosd(ref_lat)*sind(ref_long); ...
  0 cosd(ref_lat) sind(ref_lat)];

DOP_matrix_ENUP=F'*A_ECEF*F;

GDOP  = sqrt(trace(DOP_matrix_ENUP));                       % GDOP    
PDOP  = sqrt(DOP_matrix_ENUP(1,1) + DOP_matrix_ENUP(2,2) + DOP_matrix_ENUP(3,3));       % PDOP
HDOP  = sqrt(DOP_matrix_ENUP(1,1) + DOP_matrix_ENUP(2,2));                % HDOP
VDOP  = sqrt(DOP_matrix_ENUP(3,3));                         % VDOP


