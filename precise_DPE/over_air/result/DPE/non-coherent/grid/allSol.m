clear posDPEECEF
%%
% load("posDPEECEF_01_02.mat")
% load("posDPEECEF_03_04.mat")
% load("posDPEECEF_05_06.mat")
% load("posDPEECEF_07_08.mat")
% load("posDPEECEF_09_10.mat")
% load("posDPEECEF_11_12.mat")

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
pdpe_ars_20ms = load('result\PDPE\ars0321\posPDPEECEF_20ms_ars.mat');
posPDPEECEF_ars_20ms = pdpe_ars_20ms.posDPEECEF;
% pdpe_ars_20ms = load('result\PDPE\ars0325\posPDPEECEF_20ms_ars.mat');
% posPDPEECEF_ars_20ms = pdpe_ars_20ms.posPDPEECEF;


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

% mean_pos2SPECEF = mean(pos2SPECEF(1:numSol,:), 1);
% mean_posDPEECEF = mean(posDPEECEF_ars_1ms(1:numSol,:), 1);
mean_pos2SPENU = mean(pos2SPENU(1:numSol,:), 1);
mean_posDPEENU_ars_1ms = mean(posDPEENU_ars_1ms(1:numSol,:), 1);
mean_posDPEENU_ars_20ms = mean(posDPEENU_ars_20ms(1:numSol,:), 1);
mean_posPDPEENU_ars_20ms = mean(posPDPEENU_ars_20ms(1:numSol,:), 1);

error2SP=sqrt((mean_pos2SPENU(1)-settings.truePosition.E)^2 ...
    +(mean_pos2SPENU(2)-settings.truePosition.N)^2 ...
    +(mean_pos2SPENU(3)-settings.truePosition.U)^2)
errorDPE_ars_1ms=sqrt((mean_posDPEENU_ars_1ms(1)-settings.truePosition.E)^2 ...
    +(mean_posDPEENU_ars_1ms(2)-settings.truePosition.N)^2 ...
    +(mean_posDPEENU_ars_1ms(3)-settings.truePosition.U)^2)
errorDPE_ars_20ms=sqrt((mean_posDPEENU_ars_20ms(1)-settings.truePosition.E)^2 ...
    +(mean_posDPEENU_ars_20ms(2)-settings.truePosition.N)^2 ...
    +(mean_posDPEENU_ars_20ms(3)-settings.truePosition.U)^2)
errorPDPE_ars_20ms=sqrt((mean_posPDPEENU_ars_20ms(1)-settings.truePosition.E)^2 ...
    +(mean_posPDPEENU_ars_20ms(2)-settings.truePosition.N)^2 ...
    +(mean_posPDPEENU_ars_20ms(3)-settings.truePosition.U)^2)
% plot
% solPlot(truePosECEF, pos2SPECEF(1:numSol,:), posDPEECEF(1:numSol,:))

figure(111);
plot3(navSolutions.E(1:numSol) - settings.truePosition.E, ...
      navSolutions.N(1:numSol) - settings.truePosition.N, ... 
      navSolutions.U(1:numSol)- settings.truePosition.U, '+',...
      'LineWidth', 1.5, 'Color',"#0072BD");
hold on;
plot3(posDPEENU_ars_1ms(1:numSol,1) - settings.truePosition.E, ...
      posDPEENU_ars_1ms(1:numSol,2) - settings.truePosition.N, ... 
      posDPEENU_ars_1ms(1:numSol,3) - settings.truePosition.U, 'o',...
      'LineWidth', 1.5, 'Color',"#D95319");
plot3(posDPEENU_ars_20ms(1:numSol,1) - settings.truePosition.E, ...
      posDPEENU_ars_20ms(1:numSol,2) - settings.truePosition.N, ... 
      posDPEENU_ars_20ms(1:numSol,3) - settings.truePosition.U, 'o',...
      'LineWidth', 1.5, 'Color',"#77AC30");
plot3(posPDPEENU_ars_20ms(1:numSol,1) - settings.truePosition.E, ...
      posPDPEENU_ars_20ms(1:numSol,2) - settings.truePosition.N, ... 
      posPDPEENU_ars_20ms(1:numSol,3) - settings.truePosition.U, 'o',...
      'LineWidth', 1.5, 'Color',"#7205F5");
                  
%Plot the reference point
plot3(0, 0, 0, 'r+', 'LineWidth', 2.5, 'MarkerSize', 15);
hold off;

view  (0, 90);
axis  ('equal');
grid  ('minor');    

text2SP = ['2SP noncoherent 1ms (mean error = ' num2str(error2SP) ' [m])'];
textDPE_ars_1ms = ['DPE noncoherent 1ms (mean error = ' num2str(errorDPE_ars_1ms) ' [m])'];
textDPE_ars_20ms = ['DPE noncoherent 20ms (mean error = ' num2str(errorDPE_ars_20ms) ' [m])'];
textPDPE_ars_20ms = ['PDPE noncoherent 20ms (mean error = ' num2str(errorPDPE_ars_20ms) ' [m])'];

xlim([-30,30])
ylim([-30,30])
legend(text2SP, textDPE_ars_1ms, textDPE_ars_20ms, textPDPE_ars_20ms, 'True Position');
title ('Positions in UTM system');
xlabel('East (m)');
ylabel('North (m)');
zlabel('Upping (m)');
