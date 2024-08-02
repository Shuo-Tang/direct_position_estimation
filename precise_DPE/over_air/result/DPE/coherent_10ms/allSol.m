clear posDPEECEF
%%
% load("posDPEECEF_01_02.mat")
% load("posDPEECEF_03_04.mat")
% load("posDPEECEF_05_06.mat")
% load("posDPEECEF_07_08.mat")
% load("posDPEECEF_09_10.mat")
% load("posDPEECEF_11_12.mat")

%%
numFiles = 38;
sol = cell(numFiles, 1);
listing = dir();
for k = 1:numFiles
    sol{k} = load(listing(k+3).name);
end

posDPEECEF = zeros(75, 3);
for k = 1:numFiles-1
    posDPEECEF(2*k-1:2*k,:) = sol{k, 1}.posDPEECEF(2*k-1:2*k,:);
end
posDPEECEF(75,:) = sol{38, 1}.posDPEECEF(75,:);

save('posDPEECEF_01_75.mat', 'posDPEECEF')

%%
addpath include             % The software receiver functions
addpath geoFunctions        % Position calculation related functions 
settings = initSettings();
numSol = 75;
load("loadData/navSolutions_1ms.mat")
truePosECEF = [4777973.177, 176346.307, 4207663.62];
pos2SPECEF = [navSolutions.X',navSolutions.Y',navSolutions.Z'];
pos2SPENU = [navSolutions.E',navSolutions.N',navSolutions.U'];
posDPEENU = zeros(numSol, 3);
for k = 1:numSol
    [posDPEENU(k, 1), posDPEENU(k, 2),...
        posDPEENU(k, 3)] = cart2utm(posDPEECEF(k, 1), posDPEECEF(k, 2),...
        posDPEECEF(k, 3), navSolutions.utmZone);
end

mean_pos2SPECEF = mean(pos2SPECEF(1:numSol,:), 1);
mean_posDPEECEF = mean(posDPEECEF(1:numSol,:), 1);
mean_pos2SPENU = mean(pos2SPENU(1:numSol,:), 1);
mean_posDPEENU = mean(posDPEENU(1:numSol,:), 1);

error2SP=sqrt((mean_pos2SPENU(1)-settings.truePosition.E)^2 ...
    +(mean_pos2SPENU(2)-settings.truePosition.N)^2 ...
    +(mean_pos2SPENU(3)-settings.truePosition.U)^2)
errorDPE=sqrt((mean_posDPEENU(1)-settings.truePosition.E)^2 ...
    +(mean_posDPEENU(2)-settings.truePosition.N)^2 ...
    +(mean_posDPEENU(3)-settings.truePosition.U)^2)
% plot
solPlot(truePosECEF, pos2SPECEF(1:numSol,:), posDPEECEF(1:numSol,:))

figure(110);
plot3(navSolutions.E(1:numSol) - settings.truePosition.E, ...
      navSolutions.N(1:numSol) - settings.truePosition.N, ... 
      navSolutions.U(1:numSol)- settings.truePosition.U, '+',...
      'LineWidth', 1.5, 'Color',"#0072BD");
hold on;
plot3(posDPEENU(1:numSol,1) - settings.truePosition.E, ...
      posDPEENU(1:numSol,2) - settings.truePosition.N, ... 
      posDPEENU(1:numSol,3) - settings.truePosition.U, 'o',...
      'LineWidth', 1.5, 'Color',"#D95319");
                  
%Plot the reference point
plot3(0, 0, 0, 'r+', 'LineWidth', 2.5, 'MarkerSize', 15);
hold off;

view  (0, 90);
axis  ('equal');
grid  ('minor');    

text2SP = ['2SP Estimation (Mean 3D error = ' num2str(error2SP) ' [m])'];
textDPE = ['DPE Estimation (Mean 3D error = ' num2str(errorDPE) ' [m])'];

xlim([-30,30])
ylim([-30,30])
legend(text2SP, textDPE, 'True Position')%, 'DPE Estimation', 'True Position');
title ('Positions in UTM system');
xlabel('East (m)');
ylabel('North (m)');
zlabel('Upping (m)');
