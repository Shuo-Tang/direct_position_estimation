clear
clc
close all
%% CNo
% load("result/monte_carlo/pos_error_CNo_DPE.mat")
% load("result/monte_carlo/pos_error_CNo_SD.mat")
% cno_candidates = (30:1:50);
% figure(1)
% plot(cno_candidates, pos_error_CNo_DPE, "LineWidth", 2, "Color", "#77AC30")
% hold on
% plot(cno_candidates, pos_error_CNo_SD, "LineWidth", 2, "Color", "#0072BD")
% xlabel("CNo (dB-Hz)", 'Interpreter', 'latex', 'FontSize',16)
% ylabel("Positioning Error (m)", 'Interpreter', 'latex', 'FontSize',16)
% legend("DPE", "DPE-SD", 'Interpreter', 'latex', 'FontSize',16)

%% Sigma range
% load("result/monte_carlo/pos_error_sigma_DPE.mat")
% load("result/monte_carlo/pos_error_sigma_SD.mat")
% sigma_candidates = (0.1:0.5:12);
% figure(2)
% plot(sigma_candidates, pos_error_sigma_DPE, "LineWidth", 2, "Color", "#77AC30")
% hold on
% plot(sigma_candidates, pos_error_sigma_SD, "LineWidth", 2, "Color", "#0072BD")
% xlabel("Pesudorange Uncertainty $\sigma_\rho$ (m)", 'Interpreter', 'latex', 'FontSize',16)
% ylabel("Positioning Error (m)", 'Interpreter', 'latex', 'FontSize',16)
% legend("DPE", "DPE-SD", 'Interpreter', 'latex', 'FontSize',16)

%% iono error
% load("result/monte_carlo/pos_error_iono_DPE.mat")
% load("result/monte_carlo/pos_error_iono_SD.mat")
% iono_candidates = (0.5:0.2:15);
% figure(3)
% plot(iono_candidates, pos_error_iono_DPE, "LineWidth", 2, "Color", "#77AC30")
% hold on
% plot(iono_candidates, pos_error_iono_SD, "LineWidth", 2, "Color", "#0072BD")
% xlabel("Ionospheric Error (m)", 'Interpreter', 'latex', 'FontSize',16)
% ylabel("Positioning Error (m)", 'Interpreter', 'latex', 'FontSize',16)
% legend("DPE", "DPE-SD", 'Interpreter', 'latex', 'FontSize',16)

%%
load("result/real_data/sol_ecef_2sp.mat")
load("result/real_data/sol_ecef_dpe.mat");
load("result/real_data/sol_ecef_dpe_sd.mat");
X = 4777973.177; Y = 176346.307; Z = 4207663.62;
settings.truePosition.ECEF = [X, Y, Z];
LLA = ecef2lla([X, Y, Z]);
settings.truePosition.LLA = LLA; 
utmZone = 31;
[E, N, U] = cart2utm(X, Y, Z, utmZone);

nSol_2sp = size(sol_ecef_2sp, 1);
sol_enu_2sp = zeros(nSol_2sp, 3);
for iSol = 1:nSol_2sp
    [sol_enu_2sp(iSol, 1), sol_enu_2sp(iSol, 2), sol_enu_2sp(iSol, 3)] =...
        cart2utm(sol_ecef_2sp(iSol, 1), sol_ecef_2sp(iSol, 2),...
        sol_ecef_2sp(iSol, 3), utmZone);
end
nSol_dpe = size(sol_ecef_dpe, 1);
sol_enu_dpe = zeros(nSol_dpe, 3);
sol_enu_dpe_sd = zeros(nSol_dpe, 3);
for iSol = 1:nSol_dpe
    [sol_enu_dpe(iSol, 1), sol_enu_dpe(iSol, 2), sol_enu_dpe(iSol, 3)] =...
        cart2utm(sol_ecef_dpe(iSol, 1), sol_ecef_dpe(iSol, 2),...
        sol_ecef_dpe(iSol, 3), utmZone);
    [sol_enu_dpe_sd(iSol, 1), sol_enu_dpe_sd(iSol, 2), sol_enu_dpe_sd(iSol, 3)] =...
        cart2utm(sol_ecef_dpe_sd(iSol, 1), sol_ecef_dpe_sd(iSol, 2),...
        sol_ecef_dpe_sd(iSol, 3), utmZone);
end

% compute error
truePosENU = [E, N, U];
error_2sp_3d = sol_enu_2sp - repmat(truePosENU, nSol_2sp, 1);
error_2sp = vecnorm(error_2sp_3d, 2, 2);
error_dpe_3d = sol_enu_dpe - repmat(truePosENU, nSol_dpe, 1);
error_dpe = vecnorm(error_dpe_3d, 2, 2);
error_dpe_sd_3d = sol_enu_dpe_sd - repmat(truePosENU, nSol_dpe, 1);
error_dpe_sd = vecnorm(error_dpe_sd_3d, 2, 2);

mean_2sp = mean(sol_enu_2sp, 1);
mean_dpe = mean(sol_enu_dpe, 1);
mean_dpe_sd = mean(sol_enu_dpe_sd, 1);
mean_error_dpe = norm(mean_2sp - truePosENU);
mean_error_2sp = norm(mean_dpe - truePosENU);
mean_error_dpe_sd = norm(mean_dpe_sd - truePosENU);

q75th_2sp = quantile(error_2sp, 0.25);
q75th_dpe = quantile(error_dpe, 0.25);
q75th_dpe_sd = quantile(error_dpe_sd, 0.25);
[x_2sp_75th, y_2sp_75th] = circle([mean_2sp(1)-E, mean_2sp(2)-N], q75th_2sp);
[x_dpe_75th, y_dpe_75th] = circle([mean_dpe(1)-E, mean_dpe(2)-N], q75th_dpe) ;
[x_dpe_sd_75th, y_dpe_sd_75th] = circle([mean_dpe_sd(1)-E, mean_dpe_sd(2)-N], q75th_dpe_sd);



figure(112);
plot3(sol_enu_2sp(:,1) - E, ...
      sol_enu_2sp(:,2) - N, ... 
      sol_enu_2sp(:,3) - U, 'o',...
      'LineWidth', 1.5, 'Color',"#7205F5", "MarkerSize", 3);
hold on;
plot3(sol_enu_dpe(:,1) - E, ...
      sol_enu_dpe(:,2) - N, ... 
      sol_enu_dpe(:,3) - U, 'o',...
      'LineWidth', 1.5, 'Color',"#048704", "MarkerSize", 3);
plot3(sol_enu_dpe_sd(:,1) - E, ...
      sol_enu_dpe_sd(:,2) - N, ... 
      sol_enu_dpe_sd(:,3) - U, 'o',...
      'LineWidth', 1.5, 'Color',"#4DBEEE", "MarkerSize", 3);

% Plot the mean estimation and reference point
plot3(mean_2sp(1) - E, mean_2sp(2) - N, mean_2sp(3) - U, '+', ...
    'LineWidth', 2.5, 'MarkerSize', 15, 'Color',"#7205F5");
plot3(mean_dpe(1) - E, mean_dpe(2) - N, mean_dpe(3) - U, '+',...
    'LineWidth', 2.5, 'MarkerSize', 15, 'Color',"#048704");
plot3(mean_dpe_sd(1) - E, mean_dpe_sd(2) - N, mean_dpe_sd(3) - U, '+',...
    'LineWidth', 2.5, 'MarkerSize', 15, 'Color',"#4DBEEE");


% plot 75th quantile
z = zeros(1, length(x_dpe_75th));
plot3(x_2sp_75th, y_2sp_75th, z, '-.', 'LineWidth', 1.5, 'Color',"#7205F5");%"#41ff12"
plot3(x_dpe_75th, y_dpe_75th, z, '-.', 'LineWidth', 1.5, 'Color',"#048704");%"#41ff12"
plot3(x_dpe_sd_75th, y_dpe_sd_75th, z, '-.', 'LineWidth', 1.5, 'Color',"#4DBEEE");

plot3(0, 0, 0, 'r+', 'LineWidth', 2.5, 'MarkerSize', 15);
hold off;
% set view
view  (0, 90);
axis  ('equal');
grid  ('minor'); 
% legend
text_2sp = ['Mean of 2SP Estimation (Error = ' num2str(mean_error_2sp) ' [m])'];
text_dpe = ['Mean of DPE Estimation (Error = ' num2str(mean_error_dpe) ' [m])'];
text_dpe_sd = ['Mean of DPE-SD Estimation (Error = ' num2str(mean_error_dpe_sd) ' [m])'];
text_2sp_75thq = '2SP Error 75th Quantile Centered at Mean';
text_dpe_75thq = 'DPE Error 75th Quantile Centered at Mean';
text_dpe_sd_75thq = 'DPE-SD Error 75th Quantile Centered at Mean';
h_l = legend("2SP", "DPE", "DPE-SD", text_2sp, text_dpe, text_dpe_sd,...
    text_2sp_75thq, text_dpe_75thq, text_dpe_sd_75thq, 'True Position', FontSize=14);
set(h_l, 'Location', 'eastoutside');
ax = gca;
ax.Position = [0.1, 0.1, 0.7, 0.8]; % [left, bottom, width, height]
title ('Positions Estimate in UTM system', FontSize=14);
xlabel('East (m)', FontSize=14);
ylabel('North (m)', FontSize=14);
zlabel('Upping (m)');
function [x, y] = circle(center,r)
ang=0:0.01:2*pi; 
x = center(1) + r*cos(ang);
y = center(2) + r*sin(ang);
% plot(center(1)+x,center(2)+y);
end