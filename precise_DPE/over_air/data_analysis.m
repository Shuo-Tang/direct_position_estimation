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
numSol = 62;
load("loadData/navSolutions_1ms.mat")
dep_ars_1ms = load('result\DPE\non-coherent\ars\posDPEECEF_ars.mat');
posDPEECEF_ars_1ms = dep_ars_1ms.posDPEECEF;
dpe_ars_20ms = load('result\DPE\coherent_20ms\ars\posDPEECEF_20ms_ars.mat');
posDPEECEF_ars_20ms = dpe_ars_20ms.posDPEECEF;
% pdpe_ars_20ms = load('result\PDPE\ars0321\posPDPEECEF_20ms_ars.mat');
% posPDPEECEF_ars_20ms = pdpe_ars_20ms.posDPEECEF;
pdpe_ars_20ms = load('result\PDPE\ars0325\posPDPEECEF_20ms_ars.mat');
posPDPEECEF_ars_20ms = pdpe_ars_20ms.posPDPEECEF;
pdpe_ars_40ms = load('result\PDPE\ars0420_40ms\posPDPEECEF_40ms_ars.mat');
posPDPEECEF_ars_40ms = pdpe_ars_40ms.posPDPEECEF;

% load benchmark
truePosECEF = [4777973.177, 176346.307, 4207663.62];
pos2SPECEF = [navSolutions.X',navSolutions.Y',navSolutions.Z'];
pos2SPENU = [navSolutions.E',navSolutions.N',navSolutions.U'];

% convert to ENU
posDPEENU_ars_1ms = zeros(numSol, 3);
for k = 1:numSol
    [posDPEENU_ars_1ms(k, 1), posDPEENU_ars_1ms(k, 2),...
        posDPEENU_ars_1ms(k, 3)] = cart2utm(posDPEECEF_ars_1ms(k, 1), posDPEECEF_ars_1ms(k, 2),...
        posDPEECEF_ars_1ms(k, 3), navSolutions.utmZone);
end
posDPEENU_ars_20ms = zeros(numSol, 3);
for k = 1:numSol
    [posDPEENU_ars_20ms(k, 1), posDPEENU_ars_20ms(k, 2),...
        posDPEENU_ars_20ms(k, 3)] = cart2utm(posDPEECEF_ars_20ms(k, 1), posDPEECEF_ars_20ms(k, 2),...
        posDPEECEF_ars_20ms(k, 3), navSolutions.utmZone);
end
posPDPEENU_ars_20ms = zeros(numSol, 3);
for k = 1:numSol
    [posPDPEENU_ars_20ms(k, 1), posPDPEENU_ars_20ms(k, 2),...
        posPDPEENU_ars_20ms(k, 3)] = cart2utm(posPDPEECEF_ars_20ms(k, 1), posPDPEECEF_ars_20ms(k, 2),...
        posPDPEECEF_ars_20ms(k, 3), navSolutions.utmZone);
end

posPDPEENU_ars_40ms = zeros(numSol, 3);
for k = 1:numSol
    [posPDPEENU_ars_40ms(k, 1), posPDPEENU_ars_40ms(k, 2),...
        posPDPEENU_ars_40ms(k, 3)] = cart2utm(posPDPEECEF_ars_40ms(k, 1), posPDPEECEF_ars_40ms(k, 2),...
        posPDPEECEF_ars_40ms(k, 3), navSolutions.utmZone);
end

% mean_pos2SPECEF = mean(pos2SPECEF(1:numSol,:), 1);
% mean_posDPEECEF = mean(posDPEECEF_ars_1ms(1:numSol,:), 1);
mean_pos2SPENU = mean(pos2SPENU(1:numSol,:), 1);
mean_posDPEENU_ars_1ms = mean(posDPEENU_ars_1ms(1:numSol,:), 1);
mean_posDPEENU_ars_20ms = mean(posDPEENU_ars_20ms(1:numSol,:), 1);
mean_posPDPEENU_ars_20ms = mean(posPDPEENU_ars_20ms(1:numSol,:), 1);
mean_posPDPEENU_ars_40ms = mean(posPDPEENU_ars_40ms(1:numSol,:), 1);

%% Bias Error
% errorBias2SP=sqrt((mean_pos2SPENU(1)-settings.truePosition.E)^2 ...
%     +(mean_pos2SPENU(2)-settings.truePosition.N)^2 ...
%     +(mean_pos2SPENU(3)-settings.truePosition.U)^2)
% errorBiasDPE_ars_1ms=sqrt((mean_posDPEENU_ars_1ms(1)-settings.truePosition.E)^2 ...
%     +(mean_posDPEENU_ars_1ms(2)-settings.truePosition.N)^2 ...
%     +(mean_posDPEENU_ars_1ms(3)-settings.truePosition.U)^2)
errorBiasDPE_ars_20ms = sqrt((mean_posDPEENU_ars_20ms(1)-settings.truePosition.E)^2 ...
    +(mean_posDPEENU_ars_20ms(2)-settings.truePosition.N)^2 ...
    +(mean_posDPEENU_ars_20ms(3)-settings.truePosition.U)^2)
errorBiasPDPE_ars_20ms = sqrt((mean_posPDPEENU_ars_20ms(1)-settings.truePosition.E)^2 ...
    +(mean_posPDPEENU_ars_20ms(2)-settings.truePosition.N)^2 ...
    +(mean_posPDPEENU_ars_20ms(3)-settings.truePosition.U)^2)
errorBiasPDPE_ars_40ms = sqrt((mean_posPDPEENU_ars_40ms(1)-settings.truePosition.E)^2 ...
    +(mean_posPDPEENU_ars_40ms(2)-settings.truePosition.N)^2 ...
    +(mean_posPDPEENU_ars_40ms(3)-settings.truePosition.U)^2)

%% RMSE
truePosENU = [settings.truePosition.E, settings.truePosition.N, settings.truePosition.U];
error_DPE_ars_20ms_3d = posDPEENU_ars_20ms(1:numSol,:) - repmat(truePosENU, numSol, 1);
error_DPE_ars_20ms = vecnorm(error_DPE_ars_20ms_3d, 2, 2);
RMSE_DPE_ars_20ms = sqrt(sum(vecnorm(error_DPE_ars_20ms_3d, 2, 2).^2)/numSol);
error_PDPE_ars_20ms_3d = posPDPEENU_ars_20ms(1:numSol,:) - repmat(truePosENU, numSol, 1);
error_PDPE_ars_20ms = vecnorm(error_PDPE_ars_20ms_3d, 2, 2);
RMSE_PDPE_ars_20ms = sqrt(sum(vecnorm(error_PDPE_ars_20ms_3d, 2, 2).^2)/numSol);
error_PDPE_ars_40ms_3d = posPDPEENU_ars_40ms(1:numSol,:) - repmat(truePosENU, numSol, 1);
error_PDPE_ars_40ms = vecnorm(error_PDPE_ars_40ms_3d, 2, 2);
RMSE_PDPE_ars_40ms = sqrt(sum(vecnorm(error_PDPE_ars_40ms_3d, 2, 2).^2)/numSol);

%% 50th and 95th Quantiles
q50th_DPE_ars_20ms = quantile(error_DPE_ars_20ms, 0.5);
q90th_DPE_ars_20ms = quantile(error_DPE_ars_20ms, 0.1);
q50th_PDPE_ars_20ms = quantile(error_PDPE_ars_20ms, 0.5);
q90th_PDPE_ars_20ms = quantile(error_PDPE_ars_20ms, 0.1);
q50th_PDPE_ars_40ms = quantile(error_PDPE_ars_40ms, 0.5);
q90th_PDPE_ars_40ms = quantile(error_PDPE_ars_40ms, 0.1);
%% Plot
figure(111);
% plot3(navSolutions.E(1:numSol) - settings.truePosition.E, ...
%       navSolutions.N(1:numSol) - settings.truePosition.N, ... 
%       navSolutions.U(1:numSol)- settings.truePosition.U, '+',...
%       'LineWidth', 1.5, 'Color',"#0072BD");
% hold on;
% plot3(posDPEENU_ars_1ms(1:numSol,1) - settings.truePosition.E, ...
%       posDPEENU_ars_1ms(1:numSol,2) - settings.truePosition.N, ... 
%       posDPEENU_ars_1ms(1:numSol,3) - settings.truePosition.U, 'o',...
%       'LineWidth', 1.5, 'Color',"#D95319");
% Plot the positioing result
plot3(posDPEENU_ars_20ms(1:numSol,1) - settings.truePosition.E, ...
      posDPEENU_ars_20ms(1:numSol,2) - settings.truePosition.N, ... 
      posDPEENU_ars_20ms(1:numSol,3) - settings.truePosition.U, 'o',...
      'LineWidth', 1.5, 'Color',"#77AC30", "MarkerSize", 3);
hold on;
plot3(posPDPEENU_ars_20ms(1:numSol,1) - settings.truePosition.E, ...
      posPDPEENU_ars_20ms(1:numSol,2) - settings.truePosition.N, ... 
      posPDPEENU_ars_20ms(1:numSol,3) - settings.truePosition.U, 'o',...
      'LineWidth', 1.5, 'Color',"#7205F5", "MarkerSize", 3);
plot3(posPDPEENU_ars_40ms(1:numSol,1) - settings.truePosition.E, ...
      posPDPEENU_ars_40ms(1:numSol,2) - settings.truePosition.N, ... 
      posPDPEENU_ars_40ms(1:numSol,3) - settings.truePosition.U, 'o',...
      'LineWidth', 1.5, 'Color',"#4DBEEE", "MarkerSize", 3);

% Plot the quantiles
[x_DPE_50th, y_DPE_50th] = circle([0, 0], q50th_DPE_ars_20ms);
[x_DPE_90th, y_DPE_90th] = circle([0, 0], q90th_DPE_ars_20ms);
[x_PDPE_20ms_50th, y_PDPE_20ms_50th] = circle([0, 0], q50th_PDPE_ars_20ms);
[x_PDPE_20ms_90th, y_PDPE_20ms_90th] = circle([0, 0], q90th_PDPE_ars_20ms);
[x_PDPE_40ms_50th, y_PDPE_40ms_50th] = circle([0, 0], q50th_PDPE_ars_40ms);
[x_PDPE_40ms_90th, y_PDPE_40ms_90th] = circle([0, 0], q90th_PDPE_ars_40ms);

z = zeros(1, length(x_DPE_50th));
plot3(x_DPE_50th, y_DPE_50th, z, '-.', 'LineWidth', 1.5, 'Color',"#77AC30");%"#41ff12"
% plot3(x_DPE_90th, y_DPE_90th, z, '-', 'LineWidth', 1, 'Color',"#0e34c9");
plot3(x_PDPE_20ms_50th, y_PDPE_20ms_50th, z, '-.', 'LineWidth', 1.5, 'Color',"#7205F5");%"#a483fc"
% plot3(x_PDPE_90th, y_PDPE_90th, z, '-', 'LineWidth', 1, 'Color',"#0d5702");
plot3(x_PDPE_40ms_50th, y_PDPE_20ms_50th, z, '-.', 'LineWidth', 1.5, 'Color',"#4DBEEE");
% plot3(x_PDPE_90th, y_PDPE_90th, z, '-', 'LineWidth', 1, 'Color',"#0d5702");

% Plot the reference point
plot3(0, 0, 0, 'r+', 'LineWidth', 2.5, 'MarkerSize', 15);
hold off;

view  (0, 90);
axis  ('equal');
grid  ('minor');    

% text2SP = ['2SP noncoherent 1ms (mean error = ' num2str(errorBias2SP) ' [m])'];
% textDPE_ars_1ms = ['DPE noncoherent 1ms (mean error = ' num2str(errorBiasDPE_ars_1ms) ' [m])'];
% textDPE_ars_20ms = ['DPE noncoherent 20ms (bias mean error = ' num2str(errorBiasDPE_ars_20ms) ' [m])'];
% textPDPE_ars_20ms = ['PDPE noncoherent 20ms (bias mean error = ' num2str(errorBiasPDPE_ars_20ms) ' [m])'];

textDPE_ars_20ms = ['DPE (RMSE = ' num2str(RMSE_DPE_ars_20ms) ' [m])'];
textPDPE_ars_20ms = ['PDPE 20ms (RMSE = ' num2str(RMSE_PDPE_ars_20ms) ' [m])'];
textPDPE_ars_40ms = ['PDPE 40ms (RMSE = ' num2str(RMSE_PDPE_ars_40ms) ' [m])'];
textDPE_50thq = ['DPE Error 50th Quantile = ' num2str(q50th_DPE_ars_20ms) ' [m])'];
textDPE_90thq = ['DPE Error 90th Quantile = ' num2str(q90th_DPE_ars_20ms) ' [m])'];
textPDPE_20ms_50thq = ['PDPE 20ms Error 50th Quantile = ' num2str(q50th_PDPE_ars_20ms) ' [m])'];
textPDPE_20ms_90thq = ['PDPE Error 90th Quantile = ' num2str(q90th_PDPE_ars_20ms) ' [m])'];
textPDPE_40ms_50thq = ['PDPE 40ms Error 50th Quantile = ' num2str(q50th_PDPE_ars_40ms) ' [m])'];
textPDPE_40ms_90thq = ['PDPE Error 90th Quantile = ' num2str(q90th_PDPE_ars_40ms) ' [m])'];

% xlim([-30,30])
% ylim([-30,30])
% legend(text2SP, textDPE_ars_1ms, textDPE_ars_20ms, textPDPE_ars_20ms, 'True Position');
xlim([-8,8])
ylim([-8,8])
legend(textDPE_ars_20ms, textPDPE_ars_20ms, textPDPE_ars_40ms,...
    textDPE_50thq, textPDPE_20ms_50thq, textPDPE_40ms_50thq,'True Position');
title ('Positions in UTM system');
xlabel('East (m)');
ylabel('North (m)');
zlabel('Upping (m)');

% table design
BiasError = [errorBiasDPE_ars_20ms; errorBiasPDPE_ars_20ms; errorBiasPDPE_ars_40ms];
RMSE = [RMSE_DPE_ars_20ms; RMSE_PDPE_ars_20ms; RMSE_PDPE_ars_40ms];
Error50thQuantiles = [q50th_DPE_ars_20ms; q50th_PDPE_ars_20ms; q50th_PDPE_ars_40ms];
Error90thQuantiles = [q90th_DPE_ars_20ms; q90th_PDPE_ars_20ms; q90th_PDPE_ars_40ms];
errorTable = table(BiasError, RMSE, Error50thQuantiles, Error90thQuantiles);
errorTable.Properties.RowNames =["DPE", "PDPE 20ms", "PDPE 40ms"];
display(errorTable)


%%
function [x, y] = circle(center,r)
ang=0:0.01:2*pi; 
x = center(1) + r*cos(ang);
y = center(2) + r*sin(ang);
% plot(center(1)+x,center(2)+y);
end
