function [localMax, costMax] = ars2d_1ms(usrInfo, satInfo, rawSignal, caCodesTable,...
        Doppler, corr, dmax, dmin, contraction, nIter)
%GENERALIZEDPATTERNSEARCH This function is for generalized pattern search
%optimization method in 2D scenario
%   Detailed explanation goes here
d = dmax;
posInit = repmat(usrInfo(1:2), 4, 1) + [10, 10;10 -10;-10 10;-10 -10]; 
pos_trial = zeros(4, 2);
cost_trial = zeros(4, 1);
for i = 1:4
    r = generateCAFs_1ms_corr([posInit(i,:),usrInfo(3:end)], satInfo, rawSignal,...
                        caCodesTable, Doppler, corr);
    pos2d_est = zeros(nIter, 2);
    cost = zeros(nIter, 1);
    pos2d_est(1,:) = posInit(i,:);
    cost(1) = r;
    J_ant = r;
    
    %% ARS Alogorithm
    for it = 1: nIter-1
        % draw random movement
        rand_pos = pos2d_est(it,:) + d*(2*rand(2,1)-1)';
        % compute cost of the random point
        usrInfo1 = usrInfo;
        usrInfo1(1:2) = rand_pos;
        J = generateCAFs_1ms_corr(usrInfo1, satInfo, rawSignal,...
                            caCodesTable, Doppler, corr);
    
        % whether it is a better point
        if J > J_ant
            pos2d_est(it+1,:) = rand_pos;
            cost(it+1) = J;
            J_ant = J;
            d = dmax;
        else
            pos2d_est(it+1,:) = pos2d_est(it,:);
            cost(it+1) = cost(it);
            d = d/contraction;
        end
        % whether the step size is too small
        if d < dmin
            d = dmax;
        end  
    end
    pos_trial(i, :) = pos2d_est(end,:);
    cost_trial(i) = cost(end);
    %% plot path
    if ~exist('fig_path')
        fig_path = open('figs/gepsTest/abs_caf_middle_28.fig');
        hold on
    end
    truePosECEF = [4777973.177, 176346.307, 4207663.62];
    plot3(pos2d_est(:,1) - truePosECEF(1), pos2d_est(:, 2) - truePosECEF(2),...
        cost, 'LineWidth', 2)
end
[costMax, maxInd] = max(cost_trial);
localMax = pos_trial(maxInd,:);


%% plot search path
% open('figs/gepsTest/abs_caf_middle_28.fig')
% hold on
% truePosECEF = [4777973.177, 176346.307, 4207663.62];
% plot3(pos2d_est(:,1) - truePosECEF(1), pos2d_est(:, 2) - truePosECEF(2),...
%     cost, 'LineWidth', 2, 'Color', "#A2142F")
end





