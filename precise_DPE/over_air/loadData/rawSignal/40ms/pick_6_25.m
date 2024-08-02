clear
clc
close all
%%
for iSol = 1:75
    loadfile = sprintf("rawSignal_sol%d.mat", iSol);
    rawSignal_40ms = load(loadfile);
    
    rawSignal = rawSignal_40ms.rawSignal(:, 25375*(6-1) + 1: 25375*25);
    savefile = sprintf("40ms_6_25/rawSignal_sol%d.mat", iSol);
    save(savefile, "rawSignal");
end