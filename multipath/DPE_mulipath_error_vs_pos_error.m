clear
close all
clc;

%% Common settings for 2SP and DPE
addpath('geoFunctions');
addpath('include');
settings = initSettings();
nSat = 8;
prn = [10, 27 ,8, 11, 16, 23, 25, 26];
phase = zeros(1, 8);
noise_type = "real";
mulipath_error_candidates = 10%(0:10:500);
nCandidates = length(mulipath_error_candidates);
use_2SP = 0;
use_DPE = 1;
%% 2SP Tracking and Solution
if use_2SP == 1
% === Settings before tracking =====
LOSdelay = [0.0704486468067703, 0.0677269105612104, 0.0699321950136609,...
         0.0740099341247871, 0.0733128195424180, 0.0743388526652213,...
         0.0685682814518592, 0.0697029249098153];
channel = generateChannelInfo(prn, phase, LOSdelay, settings);
plot_tracking = 0;
plot_navigation = 0;
save_result = 0;

pos_error_2SP = zeros(1, nCandidates);
% === Loop for different mulipath errors =====
for iCandidates = 1: nCandidates
    multipath_error = mulipath_error_candidates(iCandidates);
    % load_file_path = sprintf("load_data/RxSignal_8sat_20_46MHz_multipath_error_%d_6sat.dat",...
    %         multipath_error);
    if multipath_error == 0
        load_file_path = sprintf("load_data/RxSignal_8sat_20_46MHz_multipath_error_0.dat");
    else
        load_file_path = sprintf("load_data/RxSignal_8sat_20_46MHz_multipath_error_%d_4sat.dat",...
            multipath_error);
    end
    % === Tracking stage =====
    settings.fileName = load_file_path;
    fprintf("Reading File %s ...\n", settings.fileName)
    [fid, message] = fopen(settings.fileName, 'rb');
    [trackResults, channel] = trackingL1(fid, channel, settings, noise_type);
    fclose(fid);
    if plot_tracking == 1
        plotTracking(1:settings.numberOfChannels, trackResults, settings);
    end
    % === Solve PVT =====
    sat_info = load('load_data/sat8_info_multipath_fix_pos.mat');
    sat_pos = sat_info.sat_info(:, 1:3, :);
    clear sat_info
    navSolutions = postNavigation(trackResults, settings, sat_pos, LOSdelay);
    if plot_navigation == 1
        plotNavigation(navSolutions, settings,1);
    end
    % === Save result and compute error =====
    if save_result == 1
        save("result\multipathErrorVSposError\trackResults_LOS.mat", "trackResults" )
        save("result\multipathErrorVSposError\navSolutions_LOS.mat", "navSolutions" )
    end
    nSol = length(navSolutions.E);
    userPos = [settings.truePosition.E,settings.truePosition.N];
    solutions = [navSolutions.E; navSolutions.N]';
    errorAll = solutions - repmat(userPos, nSol ,1);
    error_norm = vecnorm(errorAll, 2, 2);
    pos_rmse = sqrt(sum(error_norm.^2) / nSol);
    pos_error_2SP(iCandidates) = pos_rmse;

end % end for multipath errors
end % end for use_2SP

%% DPE Solution
if use_DPE == 1
% === Settings before DPE =====
% user settings
read_data = 1;
nSol = fix((settings.msToProcess - settings.skipMs) / settings.navSolPeriod);
truePosECEF = settings.truePosition.ECEF;
trueCb = 2e-8;
trueVelECEF = [0 0 0];
trueCd = 0;
trueUsrInfo = [truePosECEF, trueCb, trueVelECEF, trueCd];
samplesPerCode = settings.samplingFreq * 1e-3;
ts = 1 / settings.samplingFreq;
tc = 1 / 1.023e6; % period of one chip 
grid_search = 1;
ars_search = 0;
plot_CAF = 1;
% Grid search settings
if grid_search == 1
    % set position grid
    grid_range_pos = 40;%30;
    grid_density_pos = 1;
    origin_pos = truePosECEF;
    [xCan, yCan, zCan] = grid3d(origin_pos, grid_range_pos, grid_density_pos);
end
% ARS settings
if ars_search == 1
    % set clock bias grid
    grid_range_cb = 3e-7;%
    grid_desity_cb = 2e-8;
    origin_cb = 2e-08;
    dtCan = (origin_cb - grid_range_cb:grid_desity_cb:origin_cb + grid_range_cb - grid_desity_cb);
    nDt = length(dtCan);
    % load("dt_test.mat")
    % set z grid
    grid_range_z = 64;
    grid_desity_z = 2;
    origin_z = truePosECEF(3);
    zCan = (origin_z - grid_range_z:grid_desity_z:origin_z + grid_range_z - grid_desity_z) ;
    nz = length(zCan);
    % set ARS parameters
    dMax = 20;
    dMin = 0.01;
    contraction = 2;
    nIter = 3000;
end
% === Read satellite information =====
sat_info = load('load_data/sat8_info_multipath_fix_pos.mat');
sat_info = sat_info.sat_info(:,:,settings.skipMs + 1: end);
% === Pre-Sample C/A code =====
codeValueIndex = ceil(ts * (1:samplesPerCode) / tc);
codeValueIndex = rem(codeValueIndex, settings.codeLength);
codeValueIndex(codeValueIndex==0) = settings.codeLength;
caCodeSample = zeros(nSat, samplesPerCode); 
for iSat = 1:nSat
    caCode = generateCAcode(prn(iSat));
    caCodeSample(iSat,:) = caCode(codeValueIndex);
end

% === Loop for different mulipath errors =====
pos_error_DPE = zeros(1, nCandidates);
for iCandidates = 1: nCandidates
    % === Read raw data and satellite information =====
    multipath_error = mulipath_error_candidates(iCandidates);
    % load_file_path = sprintf("load_data/RxSignal_8sat_20_46MHz_multipath_error_%d_6sat.dat",...
    %         multipath_error);
    if multipath_error == 0
        load_file_path = sprintf("load_data/RxSignal_8sat_20_46MHz_multipath_error_0.dat");
    else
        load_file_path = sprintf("load_data/RxSignal_8sat_20_46MHz_multipath_error_%d_4sat.dat",...
            multipath_error);
    end
    settings.fileName = load_file_path;
    fprintf("Reading File %s ...\n", settings.fileName)
    if read_data == 1
        [fid, message] = fopen(settings.fileName, 'rb');
        if (settings.fileType == 8) 
            dataAdaptCoeff = 8;
        end
        fseek(fid, dataAdaptCoeff * settings.skipNumberOfBytes, 'bof');
        rawSignalDPE = readDataDPE_1ms(fid, nSol, settings);
    end
    % === Loop for each ms =====
    sol_DPE_ms = zeros(nSol, 3);
    for iSol = 1:5
        % extract raw signal and satellite information
        rawSignal = rawSignalDPE(iSol, :);
        noise = randn(1,length(rawSignal));
        rawSignal = rawSignal + 0.1*noise;
        satInfo = sat_info(:,:, (iSol - 1)*settings.navSolPeriod + 1);
        % === grid Search and Plot DPE 2D CAF (only for test)=====
        if grid_search == 1  
        % generate CAF on each grid point
        nx = length(xCan);
        ny = length(yCan);
        r = zeros(nx,ny);
        for ix = 1: nx
            for iy = 1:ny
                usrInfo = [xCan(ix), yCan(iy), trueUsrInfo(3:end)];
                r(ix, iy) = generateCAFs(usrInfo, satInfo, rawSignal, caCodeSample);
            end
        end
        rMax = max(r(:));
        [yMax,xMax] = find(r == rMax);
        pos_max = [xCan(xMax),;yCan(yMax)]';
        numMax = length(yMax); 
        randomMax = randi([1, numMax]);
        sol_ecef = [pos_max(randomMax,:), truePosECEF(3)];
        sol_lla = ecef2lla(sol_ecef);
        utmZone = findUtmZone(sol_lla(1), sol_lla(2));
        [E, N, U] = cart2utm(sol_ecef(1), sol_ecef(2), sol_ecef(3), utmZone);
        sol_DPE_ms(iSol, :) = [E, N, U]; %mean(pos_max, 1);
        if plot_CAF == 1
            % plot CAF values
            plot2DCAF(xCan, yCan, r, truePosECEF)
        end
        end
        if ars_search == 1
        % === Acceleration Random Search =====
        cost = zeros(nDt, 1);
        % --- 3d search -----
        sol_ecef_dt = zeros(nDt, 3);
        parfor idt = 1:nDt
            initial_pos = truePosECEF + 20 * randn(1, 3);%[20 * randn(1, 2), 0];
            usrInfo_initial = [initial_pos, dtCan(idt), trueUsrInfo(5:end)];
            [sol_ecef_dt(idt, :), cost(idt)] = arsSearch(usrInfo_initial, satInfo,...
                rawSignal, caCodeSample, nIter, dMax, dMin, contraction);
        end
        max_cost = max(cost);
        maxInd = find(cost == max_cost);
        max_ind = round(mean(maxInd));
        sol_ecef = sol_ecef_dt(max_ind,:);
        % --- 2d search -----
        % max_cost = 0;
        % for idt = 1%1:nDt
        %     sol_ecef_dt = zeros(nz, 3);
        %     cost_dt = zeros(nz, 1);
        %     parfor iz = 1:nz
        %         initial_pos = truePosECEF + [20 * randn(1, 2), 0];%20 * randn(1, 3);
        %         usrInfo_initial = [initial_pos(1:2), zCan(iz), dt_test(iSol)/settings.c, trueUsrInfo(5:end)];
        %         [sol_ecef_dt(iz,:), cost_dt(iz)] = arsSearch(usrInfo_initial, satInfo,...
        %             rawSignal, caCodeSample, nIter, dMax, dMin, contraction);
        %     end
        %     max_cost_dt = max(cost_dt);
        %     if max_cost_dt > max_cost
        %         maxInd = find(cost_dt == max_cost_dt);
        %         max_ind = round(mean(maxInd));
        %         sol_ecef = sol_ecef_dt(max_ind,:);
        %         max_cost = max_cost_dt;
        %     end
        % end
        sol_lla = ecef2lla(sol_ecef);
        utmZone = findUtmZone(sol_lla(1), sol_lla(2));
        [E, N, U] = cart2utm(sol_ecef(1), sol_ecef(2), sol_ecef(3), utmZone);
        sol_DPE_ms(iSol, :) = [E, N, U];
        end
    end
    userPos = [settings.truePosition.E,settings.truePosition.N];
    % userPos = [settings.truePosition.E, settings.truePosition.N, settings.truePosition.U];
    solutions = sol_DPE_ms(:,1:2);
    errorAll = solutions - repmat(userPos, nSol ,1);
    error_norm = vecnorm(errorAll, 2, 2);
    pos_rmse = sqrt(sum(error_norm.^2) / nSol);
    pos_error_DPE(iCandidates) = pos_rmse;
end % end for multipath errors
end % end for use_DPE


%% Result Analysis
% figure(1);
% plot(mulipath_error_candidates, pos_error_2SP)
% save("result/multipathErrorVSposError/pos_error_2SP_500_multipath_error.mat", "pos_error_2SP")
% save("result/multipathErrorVSposError/pos_error_DPE_500_multipath_error.mat", "pos_error_DPE")
% error_2SP_A5 = load("result/multipathErrorVSposError/pos_error_2SP_500_multipath_error.mat");
% error_2SP_A5 = error_2SP_A5.pos_error_2sp;
% error_DPE_A5 = load("result/multipathErrorVSposError/pos_error_DPE_500_multipath_error.mat");
% error_DPE_A5 = error_DPE_A5.pos_error_DPE;
% error_2SP_A3 = load("result/multipathErrorVSposError/pos_error_2SP_500_multipath_error_A3.mat");
% error_2SP_A3 = error_2SP_A3.pos_error_2SP;
% mulipath_error_candidates = (0:10:500);
% figure(1)
% plot(mulipath_error_candidates, error_2SP_A5);
% hold on
% plot(mulipath_error_candidates, error_DPE_A5)
% plot(mulipath_error_candidates, error_2SP_A3);
% 
% legend("2SP 5", "DPE 5", "2SP 3")
% nSat = 8;
% figure('Position', [100, 100, 800, 600]); % [left, bottom, width, height]
% for iSat = 1:nSat
%     correlation_path = sprintf("z%d.mat", iSat);
%     load(correlation_path)
% end
% z_sum = z1 + z2 + z3 + z4 + z5 + z6 +z7 + z8;
% shift_samples = 2000;
% plot(circshift(z_sum, shift_samples))
% y = get(gca, 'ylim'); % Get current y-axis limits
% hold on; 
% plot([shift_samples + 1, shift_samples + 1], y,...
%     'r--', 'LineWidth', 2); % Plot vertical line
% shift_samples = 2000;
% z = {z1;z2;z3;z4;z5;z6;z7;z8};
% for iSat = 1:5
%     corr = circshift(z{iSat},shift_samples);
%     if iSat == 5
%         corr = circshift(z{iSat},shift_samples - 5);
%     end
%     nSample = length(corr);
%     if nSample > 20460
%         corr = corr(1:20460);
%     end
%     if nSample < 20460
%         corr = [corr, corr(end) * ones(1, 20460 - nSample)];
%     end
%     subplot(5 , 1, iSat)
%     h1 = plot(corr);
%     y = get(gca, 'ylim');
%     hold on 
%     h2 = plot([shift_samples + 1, shift_samples + 1], y,...
%     'r--', 'LineWidth', 2);
%     xlim([shift_samples + 1 - 20, shift_samples + 1 + 20])
%     set(gca, 'XTickLabel', [], 'YTickLabel', []);
%     if iSat == 1
%         legend(h1, "CAF of Each Satellite", FontSize=16)
%     end
%     if iSat == 2
%         legend(h2, "True Delay (Aligned)", FontSize=16)
%     end
%     if iSat == 3
%         ylabel('CAF Value', 'FontSize', 16)
%     end
%     if iSat == 5
%         xlabel('Sample Index', 'FontSize', 16)
%         % xticks(shift_samples + 1);
%         % xticklabels({'True Delay'}); 
%         xticks([shift_samples + 1 - 20, shift_samples + 1, shift_samples + 1 + 20]);
%         xticklabels({'-20', '0', '20'}); 
%         set(gca, 'FontSize', 16)
%     end
% end

