function [b_est, J_ant] = arsSearch_baseline(config, rawSignal, b, ref_info, sat_info, caCode,...
        pseudorange_est, nIter, dMax, dMin, contraction)
%% Memory Allocation
est_b = zeros(nIter,3);
J_path = zeros(nIter,1);
cost = zeros(nIter,1);
%% Search Initialization
r = generateCAFs_baseline(config, rawSignal, ref_info, sat_info, caCode, pseudorange_est, b);
J_ant = r;
est_b(1,:) = b;
cost(1) = r;
J_path(1) = J_ant;
% count = 0;
%% ARS Algorithm
% compare the cost value of est_gamma and random moved gamma
d = dMax;
% hwb = waitbar(0,'ARS iterating...');
% ARS algorithm iterations
for it = 1:nIter-1       
    % waitbar(it/nIter, hwb, sprintf('ARS %d/%d completed',it,nIter))
    % compute cost of the random point
    new_b = est_b(it,:) + d*(2*rand(1,3)-1);%[d*(2*rand(1,2)-1), 0];
    r = generateCAFs_baseline(config, rawSignal, ref_info, sat_info, caCode, pseudorange_est, new_b);
    J = r;
    cost(it+1) = r;

    % select or discard the point
    if J > J_ant
        est_b(it+1,:) = new_b;
        J_ant = J;
        d = dMax;
        % count = count + 1;
    else
        est_b(it+1,:) = est_b(it,:);
        d = d/contraction;
    end
    if d < dMin
        d = dMax;
    end     
    J_path(it+1) = J_ant;
end
% close(hwb)

% DPE position estimation
b_est = est_b(end,:);
%% visualization
% figure,
% plot((1:1:Niter),position_est(:,1)-param.UserPosition(1),"r")
% hold on
% plot((1:1:Niter),zeros(1,Niter),"b")
% figure,
% plot((1:1:Niter),cb_est(:,1),"r")
% hold on
% plot((1:1:Niter),param.deltaT(1,1)*ones(1,Niter),"b")

% fig = figure;
% left_color = [0 0.4470 0.7410];
% right_color = [0.4660 0.6740 0.1880];
% set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% hold on 
% x = (1:1:Niter);
% yyaxis left
% plot(x, estPos(:,1),'Color','#0072BD','LineWidth',2)
% ylabel('x in ECEF (m)')
% yyaxis right
% plot(x, cost,'Color','#77AC30','LineWidth', 0.8)
% ylabel('Cost Function Value')
% xlabel('ARS Iteration Step')
% 
% fig = figure;
% left_color = [0 0.4470 0.7410];
% right_color = [0.4660 0.6740 0.1880];
% set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% hold on 
% x = (1:1:Niter);
% yyaxis left
% plot(x, estPos(:,1),'Color','#0072BD','LineWidth',2)
% ylabel('x in ECEF (m)')
% yyaxis right
% plot(x, curJ,'Color','#77AC30','LineWidth', 0.8)
% ylabel('Current Max Cost Value')
% xlabel('ARS Iteration Step')

%% plot path
% open("result/CAF_2d.fig")
% hold on
% true_baseline = [1000, 1000, 800];
% h_path = plot3(est_b(:,1) - true_baseline(1), est_b(:,2) - true_baseline(2), J_path,...
%     'LineWidth',3, 'Color', [0, 0, 0]);
% legend(h_path, "ARS Search Path")



