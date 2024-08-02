clear
clc
%%
X = 4777973.177; Y = 176346.307; Z = 4207663.62;
user_pos = [X, Y, Z];
load("load_data\sat6_info_original.mat")
az_end = sat_info(:, 7, end);
el_end = sat_info(:, 8, end);
prn = sat_info(:, 9, end);

figure(300);
skyPlot(az_end, el_end, prn);

extra_sat_PRN = [23; 25; 26];
extra_sat_pos_vel = sat_info(3:5, 1:6, :);
extra_sat_pos_vel_1 = mean(extra_sat_pos_vel, 1) - [0, 1e7, 1e7, 0, 0, 0];
extra_sat_pos_vel = sat_info(4:5, 1:6, :);
extra_sat_pos_vel_2 = mean(extra_sat_pos_vel, 1) + [3e6, 0, -1e7, 0, 0, 0];
extra_sat_pos_vel = sat_info(3:4, 1:6, :);
extra_sat_pos_vel_3 = mean(extra_sat_pos_vel, 1) + [0, 1e6, 1e6, 0, 0, 0];
extra_sat_pos_vel = [extra_sat_pos_vel_1;extra_sat_pos_vel_2;extra_sat_pos_vel_3];
%%
nEpoch = size(extra_sat_pos_vel, 3);
sat_info_multipath = zeros(8, 9, nEpoch);
for iEpoch = 1:nEpoch
    sat_pos = [extra_sat_pos_vel_1(:,1:3,iEpoch);
                extra_sat_pos_vel_2(:,1:3,iEpoch);
                extra_sat_pos_vel_3(:,1:3,iEpoch);];
    Az = zeros(3, 1);
    El = zeros(3, 1);
    for i = 1:3
        [Az(i), El(i), ~] = topocent(user_pos', sat_pos(i,:)' - user_pos');
    end
    extra_sat_info = [extra_sat_pos_vel(:,:,iEpoch), Az, El, extra_sat_PRN];
    sat_info_multipath(:,:,iEpoch) = [sat_info(:,:,iEpoch); extra_sat_info];
end
figure(301);
az_end = sat_info_multipath(:, 7, end);
el_end = sat_info_multipath(:, 8, end);
prn = sat_info_multipath(:, 9, end);
skyPlot(az_end, el_end, prn);

sat_info = sat_info_multipath;
% save("load_data\sat8_info_multipath.mat", "sat_info")
%%
for iEpoch = 2:nEpoch
    sat_info(:,:,iEpoch) = sat_info(:,:,1);
end
% save("load_data\sat8_info_multipath_fix_pos.mat", "sat_info")