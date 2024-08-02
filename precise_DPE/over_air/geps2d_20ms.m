function [localMax, cost] = geps2d_20ms(usrInfo, satInfo, rawSignal, caCodesTable,...
        Doppler, corr, searchD, stepSize, stepDecay, convergence)
%GENERALIZEDPATTERNSEARCH This function is for generalized pattern search
%optimization method in 2D scenario
%   Detailed explanation goes here
step0 = stepSize;
r = generateCAFs_pdpe_20ms_corr(usrInfo, satInfo, rawSignal,...
                        caCodesTable, Doppler, corr);
path = [usrInfo(1:2), r];
while stepSize > convergence
    improved  = false;
    for i = 1:size(searchD,2)
        usrInfo1 = usrInfo;
        usrInfo1(1:2) = usrInfo1(1:2) + stepSize*searchD(:,i)';
        r1 = generateCAFs_pdpe_20ms_corr(usrInfo1, satInfo, rawSignal,...
                        caCodesTable, Doppler, corr);
        if r1 > r
            usrInfo = usrInfo1;
            r = r1;
            path = [path; usrInfo(1:2), r];
            improved = true;
            stepSize = step0;
            D = searchD(:,i);
            searchD(:,i) = [];
            searchD = [D, searchD];
            break
        end
    end
    if ~improved
        stepSize = stepDecay*stepSize;
    end
end
localMax = usrInfo(1:2);
cost = r;

%% plot search path
% open('figs/gepsTest/real_caf_middle.fig')
% hold on
% truePosECEF = [4777973.177, 176346.307, 4207663.62];
% plot3(path(:,1) - truePosECEF(1), path(:, 2) - truePosECEF(2),...
%     path(:, 3), 'LineWidth', 2, 'Color', "#A2142F")
end





