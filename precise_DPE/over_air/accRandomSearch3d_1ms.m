function [posEst, J_ant] = accRandomSearch3d_1ms(...
                        usrInfo, satInfo, rawSignal, caCodesTable,...
                        Doppler, corr, searchD, stepSize, stepDecay, convergence)
%% memory allocation
nIter = 1000;
estPos = zeros(nIter,3);
curJ = zeros(nIter,1);
cost = zeros(nIter,1);
%% compute the cost function of est
[localMax, r] = generalizedPatternSearch3d_1ms(usrInfo, satInfo, rawSignal, caCodesTable,...
                        Doppler, corr, searchD, stepSize, stepDecay, convergence);
usrInfo(1:3) = localMax;
J_ant = r;
cost(1) = r;
curJ(1) = J_ant;


%% ARS Algorithm
dMax = 10;
dMin = 0.001;
contraction = 2;

% compare the cost value of est_gamma and random moved gamma
d = dMax;
estPos(1,:) = usrInfo(1:3);
% hwb = waitbar(0,'ARS iterating...');
% ARS algorithm iterations
for it = 1:nIter-1       
    % waitbar(it/nIter, hwb, sprintf('ARS %d/%d completed',it,nIter))
    % draw a random movement
    % rand_position = estPos(it,:) + d*(2*rand(3,1)-1)';
    % compute cost of the random point
    usrInfo(1:3) = estPos(it,:) + d*(2*rand(3,1)-1)';
    [localMax, r] = generalizedPatternSearch3d_1ms(usrInfo, satInfo, rawSignal, caCodesTable,...
                        Doppler, corr, searchD, stepSize, stepDecay, convergence);
    usrInfo(1:3) = localMax;
    J = r;
    cost(it+1) = r;

    % test
%     norm(rand_position - param.UserPosition) < norm(position_est(it,:)- param.UserPosition)
%     J > J_ant
%     rand_position - position_est(it,:)

    % select or discard the point
    if J > J_ant
        estPos(it+1,:) = usrInfo(1:3);
        J_ant = J;
        d = dMax;
    else
        estPos(it+1,:) = estPos(it,:);
        d = d/contraction;
    end
    if d < dMin
        d = dMax;
    end     
    curJ(it+1) = J_ant;
end
% close(hwb)

% DPE position estimation
posEst = estPos(it+1,:);
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
% open("figs\gepsTest\real_caf_middle.fig")
% hold on
% truePosECEF = [4777973.177, 176346.307, 4207663.62];
% plot3(estPos(:,1) - truePosECEF(1), estPos(:,2) - truePosECEF(2), curJ,'LineWidth',2)


