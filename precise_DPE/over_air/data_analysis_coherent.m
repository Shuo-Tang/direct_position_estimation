clear posDPEECEF
%%
%% integrate DPE result
% numFiles = 38;
% sol = cell(numFiles, 1);
% listing = dir();
% for k = 1:numFiles
%     sol{k} = load(listing(k+3).name);
% end
% 
% posDPEECEF = zeros(75, 3);
% for k = 1:numFiles-1
%     posDPEECEF(2*k-1:2*k,:) = sol{k, 1}.posDPEECEF(2*k-1:2*k,:);
% end
% posDPEECEF(75,:) = sol{38, 1}.posDPEECEF(75,:);
% 
% save('posDPEECEF_01_75.mat', 'posDPEECEF')
% DPE_result = load('result\DPE\coherent_20ms\ars\posDPEECEF_20ms_ars.mat');

%% integrate PDPE result
% numFiles = 15;
% sol = cell(numFiles, 1);
% listing = dir();
% for k = 1:numFiles
%     sol{k} = load(listing(k+2).name);
% end
% 
% posPDPEECEF = zeros(75, 3);
% for k = 1:numFiles
%     posPDPEECEF(5*(k-1)+1:5*k,:) = sol{k, 1}.posDPEECEF(5*(k-1)+1:5*k,:);
% end


%%
% load data
addpath include             % The software receiver functions
addpath geoFunctions        % Position calculation related functions 
settings = initSettings();
numSol = 60;
numSol_TWO = 42;
% 2SP 1ms
load("loadData/navSolutions_1ms.mat")
% 2SP 20ms coherent
TWO_coherent_20ms_LLA = readmatrix('result/2SP/posLLA_gnss_sdr_20ms_coherent.txt');
TWO_coherent_20ms_LLA = [TWO_coherent_20ms_LLA(:,2),...
    TWO_coherent_20ms_LLA(:,1), TWO_coherent_20ms_LLA(:,3)];
TWO_coherent_20ms = lla2ecef(TWO_coherent_20ms_LLA);
% DPE 1ms coherent
dep_1ms = load('result\DPE\non-coherent\ars\posDPEECEF_ars.mat');
posDPEECEF_1ms = dep_1ms.posDPEECEF;
% DPE 20ms coherent
dpe_coherent_20ms = load('result\DPE\coherent_20ms\ars\posDPEECEF_20ms_ars.mat');
posDPEECEF_coherent_20ms = dpe_coherent_20ms.posDPEECEF;
% DPE 20ms noncoherent
pdpe_noncoherent_20ms = load('result\PDPE\ars0325\posPDPEECEF_20ms_ars.mat');
posPDPEECEF_noncoherent_20ms = pdpe_noncoherent_20ms.posPDPEECEF;
% PDPE 20ms coherent
pdpe_coherent_20ms = load('result\PDPE\ars0420_40ms\posPDPEECEF_40ms_ars.mat');
posPDPEECEF_coherent_20ms = pdpe_coherent_20ms.posPDPEECEF;

% load benchmark
truePosECEF = [4777973.177, 176346.307, 4207663.62];
pos2SPECEF = [navSolutions.X',navSolutions.Y',navSolutions.Z'];
pos2SPENU = [navSolutions.E',navSolutions.N',navSolutions.U'];

% convert to ENU
posTWOENU_coherent_20ms = zeros(numSol_TWO, 3);
for k = 1:numSol_TWO
    [posTWOENU_coherent_20ms(k, 1), posTWOENU_coherent_20ms(k, 2),...
        posTWOENU_coherent_20ms(k, 3)] = cart2utm(TWO_coherent_20ms(k, 1), TWO_coherent_20ms(k, 2),...
        TWO_coherent_20ms(k, 3), navSolutions.utmZone);
end

posDPEENU_1ms = zeros(numSol, 3);
for k = 1:numSol
    [posDPEENU_1ms(k, 1), posDPEENU_1ms(k, 2),...
        posDPEENU_1ms(k, 3)] = cart2utm(posDPEECEF_1ms(k, 1), posDPEECEF_1ms(k, 2),...
        posDPEECEF_1ms(k, 3), navSolutions.utmZone);
end
posDPEENU_coherent_20ms = zeros(numSol, 3);
for k = 1:numSol
    [posDPEENU_coherent_20ms(k, 1), posDPEENU_coherent_20ms(k, 2),...
        posDPEENU_coherent_20ms(k, 3)] = cart2utm(posDPEECEF_coherent_20ms(k, 1), posDPEECEF_coherent_20ms(k, 2),...
        posDPEECEF_coherent_20ms(k, 3), navSolutions.utmZone);
end
posPDPEENU_noncoherent_20ms = zeros(numSol, 3);
for k = 1:numSol
    [posPDPEENU_noncoherent_20ms(k, 1), posPDPEENU_noncoherent_20ms(k, 2),...
        posPDPEENU_noncoherent_20ms(k, 3)] = cart2utm(posPDPEECEF_noncoherent_20ms(k, 1), posPDPEECEF_noncoherent_20ms(k, 2),...
        posPDPEECEF_noncoherent_20ms(k, 3), navSolutions.utmZone);
end

posPDPEENU_coherent_20ms = zeros(numSol, 3);
for k = 1:numSol
    [posPDPEENU_coherent_20ms(k, 1), posPDPEENU_coherent_20ms(k, 2),...
        posPDPEENU_coherent_20ms(k, 3)] = cart2utm(posPDPEECEF_coherent_20ms(k, 1), posPDPEECEF_coherent_20ms(k, 2),...
        posPDPEECEF_coherent_20ms(k, 3), navSolutions.utmZone);
end

% mean_pos2SPECEF = mean(pos2SPECEF(1:numSol,:), 1);
% mean_posDPEECEF = mean(posDPEECEF_1ms(1:numSol,:), 1);
mean_posTWOENU_coherent_20ms = mean(posTWOENU_coherent_20ms(1:numSol_TWO,:), 1);
mean_pos2SPENU = mean(pos2SPENU(1:numSol,:), 1);
mean_posDPEENU_1ms = mean(posDPEENU_1ms(1:numSol,:), 1);
mean_posDPEENU_coherent_20ms = mean(posDPEENU_coherent_20ms(1:numSol,:), 1);
mean_posPDPEENU_noncoherent_20ms = mean(posPDPEENU_noncoherent_20ms(1:numSol,:), 1);
mean_posPDPEENU_coherent_20ms = mean(posPDPEENU_coherent_20ms(1:numSol,:), 1);

%% Bias Error
% errorBias2SP=sqrt((mean_pos2SPENU(1)-settings.truePosition.E)^2 ...
%     +(mean_pos2SPENU(2)-settings.truePosition.N)^2 ...
%     +(mean_pos2SPENU(3)-settings.truePosition.U)^2)
% errorBiasDPE_1ms=sqrt((mean_posDPEENU_1ms(1)-settings.truePosition.E)^2 ...
%     +(mean_posDPEENU_1ms(2)-settings.truePosition.N)^2 ...
%     +(mean_posDPEENU_1ms(3)-settings.truePosition.U)^2)
errorBiasTWO_coherent_20ms = sqrt((mean_posTWOENU_coherent_20ms(1)-settings.truePosition.E)^2 ...
    +(mean_posTWOENU_coherent_20ms(2)-settings.truePosition.N)^2 ...
    +(mean_posTWOENU_coherent_20ms(3)-settings.truePosition.U)^2)
errorBiasDPE_coherent_20ms = sqrt((mean_posDPEENU_coherent_20ms(1)-settings.truePosition.E)^2 ...
    +(mean_posDPEENU_coherent_20ms(2)-settings.truePosition.N)^2 ...
    +(mean_posDPEENU_coherent_20ms(3)-settings.truePosition.U)^2)
errorBiasPDPE_noncoherent_20ms = sqrt((mean_posPDPEENU_noncoherent_20ms(1)-settings.truePosition.E)^2 ...
    +(mean_posPDPEENU_noncoherent_20ms(2)-settings.truePosition.N)^2 ...
    +(mean_posPDPEENU_noncoherent_20ms(3)-settings.truePosition.U)^2)
errorBiasPDPE_coherent_20ms = sqrt((mean_posPDPEENU_coherent_20ms(1)-settings.truePosition.E)^2 ...
    +(mean_posPDPEENU_coherent_20ms(2)-settings.truePosition.N)^2 ...
    +(mean_posPDPEENU_coherent_20ms(3)-settings.truePosition.U)^2)

%% RMSE
truePosENU = [settings.truePosition.E, settings.truePosition.N, settings.truePosition.U];
error_TWO_coherent_20ms_3d = posTWOENU_coherent_20ms(1:numSol_TWO,:) - repmat(truePosENU, numSol_TWO, 1);
error_TWO_coherent_20ms = vecnorm(error_TWO_coherent_20ms_3d, 2, 2);
RMSE_TWO_coherent_20ms = sqrt(sum(vecnorm(error_TWO_coherent_20ms_3d, 2, 2).^2)/numSol_TWO);
error_DPE_coherent_20ms_3d = posDPEENU_coherent_20ms(1:numSol,:) - repmat(truePosENU, numSol, 1);
error_DPE_coherent_20ms = vecnorm(error_DPE_coherent_20ms_3d, 2, 2);
RMSE_DPE_coherent_20ms = sqrt(sum(vecnorm(error_DPE_coherent_20ms_3d, 2, 2).^2)/numSol);
error_PDPE_noncoherent_20ms_3d = posPDPEENU_noncoherent_20ms(1:numSol,:) - repmat(truePosENU, numSol, 1);
error_PDPE_noncoherent_20ms = vecnorm(error_PDPE_noncoherent_20ms_3d, 2, 2);
RMSE_PDPE_noncoherent_20ms = sqrt(sum(vecnorm(error_PDPE_noncoherent_20ms_3d, 2, 2).^2)/numSol);
error_PDPE_coherent_20ms_3d = posPDPEENU_coherent_20ms(1:numSol,:) - repmat(truePosENU, numSol, 1);
error_PDPE_coherent_20ms = vecnorm(error_PDPE_coherent_20ms_3d, 2, 2);
RMSE_PDPE_coherent_20ms = sqrt(sum(vecnorm(error_PDPE_coherent_20ms_3d, 2, 2).^2)/numSol);

%% 50th and 95th Quantiles
q50th_TWO_coherent_20ms = quantile(error_TWO_coherent_20ms, 0.5);
q90th_TWO_coherent_20ms = quantile(error_TWO_coherent_20ms, 0.1);
q50th_DPE_coherent_20ms = quantile(error_DPE_coherent_20ms, 0.5);
q90th_DPE_coherent_20ms = quantile(error_DPE_coherent_20ms, 0.1);
q50th_PDPE_noncoherent_20ms = quantile(error_PDPE_noncoherent_20ms, 0.5);
q90th_PDPE_noncoherent_20ms = quantile(error_PDPE_noncoherent_20ms, 0.1);
q50th_PDPE_coherent_20ms = quantile(error_PDPE_coherent_20ms, 0.5);
q90th_PDPE_coherent_20ms = quantile(error_PDPE_coherent_20ms, 0.1);
%% Plot
figure(112);
% plot3(navSolutions.E(1:numSol) - settings.truePosition.E, ...
%       navSolutions.N(1:numSol) - settings.truePosition.N, ... 
%       navSolutions.U(1:numSol)- settings.truePosition.U, '+',...
%       'LineWidth', 1.5, 'Color',"#0072BD");
% hold on;
% plot3(posDPEENU_1ms(1:numSol,1) - settings.truePosition.E, ...
%       posDPEENU_1ms(1:numSol,2) - settings.truePosition.N, ... 
%       posDPEENU_1ms(1:numSol,3) - settings.truePosition.U, 'o',...
%       'LineWidth', 1.5, 'Color',"#D95319");
% Plot the positioing result
plot3(posTWOENU_coherent_20ms(1:numSol_TWO,1) - settings.truePosition.E, ...
      posTWOENU_coherent_20ms(1:numSol_TWO,2) - settings.truePosition.N, ... 
      posTWOENU_coherent_20ms(1:numSol_TWO,3) - settings.truePosition.U, 'o',...
      'LineWidth', 1.5, 'Color',"#7205F5", "MarkerSize", 3);
hold on;
plot3(posDPEENU_coherent_20ms(1:numSol,1) - settings.truePosition.E, ...
      posDPEENU_coherent_20ms(1:numSol,2) - settings.truePosition.N, ... 
      posDPEENU_coherent_20ms(1:numSol,3) - settings.truePosition.U, 'o',...
      'LineWidth', 1.5, 'Color',"#048704", "MarkerSize", 3);
% plot3(posPDPEENU_noncoherent_20ms(1:numSol,1) - settings.truePosition.E, ...
%       posPDPEENU_noncoherent_20ms(1:numSol,2) - settings.truePosition.N, ... 
%       posPDPEENU_noncoherent_20ms(1:numSol,3) - settings.truePosition.U, 'o',...
%       'LineWidth', 1.5, 'Color',"#7205F5", "MarkerSize", 3);
plot3(posPDPEENU_coherent_20ms(1:numSol,1) - settings.truePosition.E, ...
      posPDPEENU_coherent_20ms(1:numSol,2) - settings.truePosition.N, ... 
      posPDPEENU_coherent_20ms(1:numSol,3) - settings.truePosition.U, 'o',...
      'LineWidth', 1.5, 'Color',"#4DBEEE", "MarkerSize", 3);

% Plot the quantiles
[x_TWO_50th, y_TWO_50th] = circle([0, 0], q50th_TWO_coherent_20ms);
[x_DPE_50th, y_DPE_50th] = circle([0, 0], q50th_DPE_coherent_20ms);
[x_DPE_90th, y_DPE_90th] = circle([0, 0], q90th_DPE_coherent_20ms);
[x_PDPE_noncoherent_20ms_50th, y_PDPE_noncoherent_20ms_50th] = circle([0, 0], q50th_PDPE_noncoherent_20ms);
[x_PDPE_noncoherent_20ms_90th, y_PDPE_noncoherent_20ms_90th] = circle([0, 0], q90th_PDPE_noncoherent_20ms);
[x_PDPE_coherent_20ms_50th, y_PDPE_coherent_20ms_50th] = circle([0, 0], q50th_PDPE_coherent_20ms);
[x_PDPE_coherent_20ms_90th, y_PDPE_coherent_20ms_90th] = circle([0, 0], q90th_PDPE_coherent_20ms);

z = zeros(1, length(x_DPE_50th));
plot3(x_TWO_50th, y_TWO_50th, z, '-.', 'LineWidth', 1.5, 'Color',"#7205F5");%"#41ff12"
plot3(x_DPE_50th, y_DPE_50th, z, '-.', 'LineWidth', 1.5, 'Color',"#048704");%"#41ff12"
% plot3(x_DPE_90th, y_DPE_90th, z, '-', 'LineWidth', 1, 'Color',"#0e34c9");
% plot3(x_PDPE_noncoherent_20ms_50th, y_PDPE_noncoherent_20ms_50th, z, '-.', 'LineWidth', 1.5, 'Color',"#7205F5");%"#a483fc"
% plot3(x_PDPE_90th, y_PDPE_90th, z, '-', 'LineWidth', 1, 'Color',"#0d5702");
plot3(x_PDPE_coherent_20ms_50th, y_PDPE_coherent_20ms_50th, z, '-.', 'LineWidth', 1.5, 'Color',"#4DBEEE");
% plot3(x_PDPE_90th, y_PDPE_90th, z, '-', 'LineWidth', 1, 'Color',"#0d5702");

% Plot the reference point
plot3(0, 0, 0, 'r+', 'LineWidth', 2.5, 'MarkerSize', 15);
hold off;

view  (0, 90);
axis  ('equal');
grid  ('minor');    

% text2SP = ['2SP noncoherent 1ms (mean error = ' num2str(errorBias2SP) ' [m])'];
% textDPE_1ms = ['DPE noncoherent 1ms (mean error = ' num2str(errorBiasDPE_1ms) ' [m])'];
% textDPE_coherent_20ms = ['DPE coherent coherent 20ms (bias mean error = ' num2str(errorBiasDPE_coherent_20ms) ' [m])'];
% textPDPE_noncoherent_20ms = ['PDPE noncoherent 20ms (bias mean error = ' num2str(errorBiasPDPE_noncoherent_20ms) ' [m])'];

textTWO_coherent_20ms = ['2SP (RMSE = ' num2str(RMSE_TWO_coherent_20ms) ' [m])'];
textDPE_coherent_20ms = ['DPE (RMSE = ' num2str(RMSE_DPE_coherent_20ms) ' [m])'];
textPDPE_noncoherent_20ms = ['PDPE noncoherent 20ms (RMSE = ' num2str(RMSE_PDPE_noncoherent_20ms) ' [m])'];
textPDPE_coherent_20ms = ['PDPE (RMSE = ' num2str(RMSE_PDPE_coherent_20ms) ' [m])'];
textTWO_coherent_20ms_50thq = ['2SP Error 50th Quantile = ' num2str(q50th_TWO_coherent_20ms) ' [m])'];
textTWO_coherent_20ms_90thq = ['2SP Error 90th Quantile = ' num2str(q90th_TWO_coherent_20ms) ' [m])'];
textDPE_coherent_20ms_50thq = ['DPE Error 50th Quantile = ' num2str(q50th_DPE_coherent_20ms) ' [m])'];
textDPE_coherent_20ms_90thq = ['DPE Error 90th Quantile = ' num2str(q90th_DPE_coherent_20ms) ' [m])'];
textPDPE_noncoherent_20ms_50thq = ['PDPE noncoherent 20ms Error 50th Quantile = ' num2str(q50th_PDPE_noncoherent_20ms) ' [m])'];
textPDPE_noncoherent_20ms_90thq = ['PDPE noncoherent 20ms Error 90th Quantile = ' num2str(q90th_PDPE_noncoherent_20ms) ' [m])'];
textPDPE_coherent_20ms_50thq = ['PDPE Error 50th Quantile = ' num2str(q50th_PDPE_coherent_20ms,'%.4f') ' [m])'];
textPDPE_coherent_20ms_90thq = ['PDPE Error 90th Quantile = ' num2str(q90th_PDPE_coherent_20ms) ' [m])'];

% xlim([-30,30])
% ylim([-30,30])
% legend(text2SP, textDPE_1ms, textDPE_coherent_20ms, textPDPE_noncoherent_20ms, 'True Position');
xlim([-8,8])
ylim([-10,6])
legend(textTWO_coherent_20ms, textDPE_coherent_20ms, textPDPE_coherent_20ms,...
    textTWO_coherent_20ms_50thq, textDPE_coherent_20ms_50thq, textPDPE_coherent_20ms_50thq,'True Position',...
    fontsize=14);
title ('Positions Estimate in UTM system', FontSize=14);
xlabel('East (m)', FontSize=14);
ylabel('North (m)', FontSize=14);
zlabel('Upping (m)');

% table design
BiasError = [errorBiasTWO_coherent_20ms;errorBiasDPE_coherent_20ms; errorBiasPDPE_coherent_20ms];
RMSE = [RMSE_TWO_coherent_20ms;RMSE_DPE_coherent_20ms; RMSE_PDPE_coherent_20ms];
Error50thQuantiles = [q50th_TWO_coherent_20ms;q50th_DPE_coherent_20ms; q50th_PDPE_coherent_20ms];
Error90thQuantiles = [q90th_TWO_coherent_20ms;q90th_DPE_coherent_20ms; q90th_PDPE_coherent_20ms];
errorTable = table(BiasError, RMSE, Error50thQuantiles, Error90thQuantiles);
errorTable.Properties.RowNames =["2SP","DPE", "PDPE"];
display(errorTable)


%%
function [x, y] = circle(center,r)
ang=0:0.01:2*pi; 
x = center(1) + r*cos(ang);
y = center(2) + r*sin(ang);
% plot(center(1)+x,center(2)+y);
end
