clear
clc
close all

%% Set parameters and positions
config = 'ConfigFile';
eval(config)
nSamples_1ms = fs * 1e-3;
plot_delays = 1;
noise_on = 1;
iono_tropo_on = 0;
range_error_on = 0;
plot_CAFs = 1;
DPE_SD_enable = 0;
DPE_standard_enable = 1;
% === Reference station =====
ref_info = zeros(1, 8);
ref_info(1:3) = [4777973.177,176346.307,4207663.620]; % position
ref_info(4) = 0; % clock bias
ref_info(5:7) = [0, 0, 0]; % velocity
ref_info(8) = 0; % clock drift

% === Rover =====
% ref_llh = ecef2lla(ref_info(1:3));
% wgs84 = wgs84Ellipsoid('meter');
% [ref_n,ref_e,ref_d] = ecef2ned(ref_info(1),ref_info(2), ref_info(3),...
%     ref_llh(1),ref_llh(2),ref_llh(3),wgs84);
% ref_ned = [ref_n, ref_e, ref_d];
% baseline = [1000, 1000, 50];%
% rov_ned = ref_ned + baseline;
% [rov_x,rov_y,rov_z] = ned2ecef(rov_ned (1),rov_ned (2),rov_ned(3),...
%     ref_llh(1),ref_llh(2),ref_llh(3),wgs84);
% rov_ecef = [rov_x,rov_y,rov_z];
baseline = [1000, 1000, 800];%[1000, 1000, 800]
rov_ecef = ref_info(1:3) + baseline;
rov_info = zeros(1, 8);
rov_info(1:3) = rov_ecef; % position
rov_info(4) = 0; % clock bias
rov_info(5:7) = [0, 0, 0]; % velocity
rov_info(8) = 0; % clock drift

% === Satellite =====
sat_info_original = load("load_data/sat8_info_multipath_fix_pos.mat");
sat_info_original = sat_info_original.sat_info(:,:,1);
nSat = size(sat_info_original, 1);
sat_cb = zeros(nSat, 1);
sat_cd = zeros(nSat, 1);
sat_info = [sat_info_original(:,1:3), sat_cb,...
    sat_info_original(:,4:6), sat_cd, sat_info_original(:,9)];
caCode = zeros(32, 1023);
for i = 1:32
    caCode(i,:) = genCAcode(i);
end
%% Generate two signals and observations
iono_error = 10;
tropo_error = 3;
iono_sigma = 1;
tropo_sigma = 0.1;
CNo = 50;
% === Reference station =====
if iono_tropo_on == 1
    % iono_ref = iono_error + iono_sigma*randn(nSat, 1);
    % tropo_ref = tropo_error + tropo_sigma*randn(nSat, 1);
    iono_ref = load("load_data/iono_ref.mat");
    iono_ref = iono_ref.iono_ref;
    tropo_ref = load("load_data/tropo_ref.mat");
    tropo_ref = tropo_ref.tropo_ref;
else
    iono_ref = zeros(nSat, 1);
    tropo_ref = zeros(nSat, 1);
end
[x_ref, pseudorange_true_ref] =...
    signalGenerate_iono(config, ref_info, sat_info, caCode, iono_ref, tropo_ref);
if noise_on == 1
    raw_signal_ref = receivedSignal(x_ref, config, CNo);
else
    raw_signal_ref = x_ref;
end
% === Rover station ===== 
if iono_tropo_on == 1
    % iono_rov = iono_error + iono_sigma*randn(nSat, 1);
    % tropo_rov = tropo_error + tropo_sigma*randn(nSat, 1);
    iono_rov = load("load_data/iono_rov.mat");
    iono_rov = iono_rov.iono_rov;
    tropo_rov = load("load_data/tropo_rov.mat");
    tropo_rov = tropo_rov.tropo_rov;
else
    iono_rov = zeros(nSat, 1);
    tropo_rov = zeros(nSat, 1);
end
[x_rov, pseudorange_true_rov] = signalGenerate_iono(config, rov_info, sat_info, caCode, iono_rov, tropo_rov);
if noise_on == 1
    raw_signal_rov = receivedSignal(x_rov, config, CNo);
else
    raw_signal_rov = x_rov;
end
%% Correlate between two raw signals
corr_no_noise = abs(ifft(fft(x_rov).*conj(fft(x_ref))).^2);
corr_raw = abs(ifft(fft(raw_signal_rov).*conj(fft(raw_signal_ref))).^2);
delay_vector = (pseudorange_true_rov - pseudorange_true_ref) / c;
delay_ind = round(delay_vector * fs);
if plot_delays == 1    
    figure(1)
    plot_shift = 1e4;
    % plot(circshift(corr_no_noise / max(corr_no_noise), plot_shift), 'LineWidth', 2)
    % hold on
    plot(circshift(corr_raw / max(corr_raw), plot_shift), 'LineWidth', 2) % "Color", "#77AC30"
    y = ylim;
    for iSat = 1:nSat
        line([delay_ind(iSat)+plot_shift  delay_ind(iSat)+plot_shift],...
            y, 'Color', 'r', 'LineWidth', 2, 'LineStyle', ':'); % Plot vertical line
    end
    xlim([plot_shift-400, plot_shift+400])
    xticks([plot_shift-400, plot_shift+400]);
    xticklabels({'-400', '400'})
    yticks([0, 1])
    ax = gca;
    ax.XAxis.FontSize = 14; % Font size
    ax.YAxis.FontSize = 14; % Font size
    xlabel("Sample Index", 'Interpreter', 'latex', 'FontSize',16)
    ylabel("Normalized CAFs", 'Interpreter', 'latex', 'FontSize',16)
    legend('Correlation between Raw Signals', 'True delays', 'Interpreter', 'latex', 'FontSize',16)
end

%% Correlate between raw signal and reconstructed signal
if range_error_on == 1
    pseudorange_est_ref = pseudorange_true_ref + 1*randn(nSat, 1)
else
    pseudorange_est_ref = pseudorange_true_ref
end
x_ref_constructed = signalGenerate_given_each_range(config,...
                ref_info, sat_info, caCode, pseudorange_est_ref);
% x_ref_constructed_2 = signalGenerate_given_baseline(config,...
%                 ref_info, sat_info, caCode, pseudorange_est_ref, b);
% x_re
corr_constructed = zeros(nSat, nSamples_1ms);
for iSat = 1:nSat
%     x_ref_constructed(iSat,:) = signalGenerate_given_each_range(config,...
%         ref_info, sat_info, caCode, pseudorange_est_ref(iSat));
    corr_constructed(iSat,: ) = abs(ifft(fft(raw_signal_rov).*...
        conj(fft(x_ref_constructed(iSat,:)))).^2);
end
corr_constructed_sum = sum(corr_constructed, 1);
delay_vector = (pseudorange_true_rov - pseudorange_true_ref) / c;
delay_ind = round(delay_vector * fs);
if plot_delays == 1
    figure(2)
    set(gcf, 'Position', [400, 400, 550, 450]);
    for iSat = 1:6
        if iSat == 5
            continue;
        end
        subplot(6,1,iSat)
        plot_shift = 1e4;
        plot(circshift(corr_constructed(iSat,: ) / max(corr_constructed(iSat,: )),...
            plot_shift), 'LineWidth', 2, "Color", "#77AC30")
        hold on
        y = ylim;
        line([delay_ind(iSat)+plot_shift  delay_ind(iSat)+plot_shift],...
            y, 'Color', 'r', 'LineWidth', 2, 'LineStyle', ':'); % Plot vertical line
        xlim([plot_shift-400, plot_shift+400])
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        if iSat == 6
            xlabel("Sample Index", 'Interpreter', 'latex', 'FontSize',16)
        end
        if iSat == 4
            ylabel("Normalized CAF", 'Interpreter', 'latex', 'FontSize',16)
            yticks([0,1])
        end
    end     
    legend('Correlation with Reconstructed Signal', 'True delays', 'Interpreter', 'latex', 'FontSize',14)
end

%% implement DPE-SD -- grid search
if DPE_SD_enable == 1
% === settings =====
sat_pos = sat_info(:,1:3);
ref_pos = ref_info(1:3);
pos_diff = sat_pos - repmat(ref_pos, nSat, 1);
unit_direction = pos_diff ./ repmat(vecnorm(pos_diff, 2, 2), 1, 3);

search_range = 50;
search_density = 1;
[xCan, yCan] = grid2d(baseline(1:2), search_range, search_density);
nx = length(xCan);
ny = length(yCan);

% === implement DPE-SD =====
r = zeros(nx, ny);
for ix = 1:nx
    parfor iy = 1:ny
        cost = 0;
        b = [xCan(ix), yCan(iy), baseline(3)];
        x_ref_constructed = signalGenerate_given_baseline(config,...
                ref_info, sat_info, caCode, pseudorange_est_ref, b);
        for iSat = 1:nSat
            cost = cost + abs(sum(x_ref_constructed(iSat,:).*conj(raw_signal_rov)))^2;           
        end
        r(ix, iy) = cost;
    end
end
% === look for maximum =====
rMax = max(r(:));
[yMax,xMax] = find(r == rMax)
b_x = mean(xCan(xMax));
b_y = mean(yCan(yMax));
norm([b_x, b_y, baseline(3)] - baseline)
plot2DCAF(xCan, yCan, r, baseline);
end

%% implement standard DPE -- grid search
if DPE_standard_enable == 1
search_range0 = 80;
search_density0 = 2;
[xCan0, yCan0] = grid2d(rov_info(1:2), search_range0, search_density0);
nx0 = length(xCan0);
ny0 = length(yCan0);

% === implement standard DPE =====
r0 = zeros(nx0, ny0);
for ix = 1:nx0
    for iy = 1:ny0
        rov_info_dpe = [xCan0(ix), yCan0(iy), rov_info(3:end)];
        r0(ix, iy) = generateCAFs(config, raw_signal_rov, rov_info_dpe, sat_info, caCode);
    end
end
% === look for maximum =====
rMax0 = max(r0(:));
[yMax0,xMax0] = find(r0 == rMax0)
pos_x = mean(xCan0(xMax0));
pos_y = mean(yCan0(yMax0));
norm([pos_x, pos_y, rov_info(3)] - rov_info(1:3))
plot2DCAF(xCan0, yCan0, r0, rov_info(1:2));


end