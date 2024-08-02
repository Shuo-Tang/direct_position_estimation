clear
clc
close all

%% multipath error envelope
% figure(1)
% error_2SP_A5 = load("result/multipathErrorVSposError/pos_error_2SP_500_multipath_error.mat");
% error_2SP_A5 = error_2SP_A5.pos_error_2SP;
% error_DPE_A5 = load("result/multipathErrorVSposError/pos_error_DPE_500_multipath_error.mat");
% error_DPE_A5 = error_DPE_A5.pos_error_DPE;
% error_2SP_A3 = load("result/multipathErrorVSposError/pos_error_2SP_500_multipath_error_A3.mat");
% error_2SP_A3 = error_2SP_A3.pos_error_2SP;
% error_DPE_A3 = load("result/multipathErrorVSposError/pos_error_DPE_500_multipath_error_A3.mat");
% error_DPE_A3 = error_DPE_A3.pos_error_DPE;
% mulipath_error_candidates = (0:10:500);
% plot(mulipath_error_candidates, error_2SP_A5, 'Color','#0072bd', 'LineWidth', 2, 'LineStyle',':');
% hold on
% plot(mulipath_error_candidates, error_2SP_A3, 'Color','#6bcffa', 'LineWidth', 2, 'LineStyle',':');
% plot(mulipath_error_candidates, error_DPE_A5, 'Color','#c90428', 'LineWidth', 2)
% plot(mulipath_error_candidates, error_DPE_A3, 'Color','#ff6929', 'LineWidth', 2);
% 
% legend("2SP \,\,a=0.5", "2SP \,\,a=0.3", "DPE a=0.5", "DPE a=0.3",...
%     'Interpreter', 'latex', 'FontSize',14)
% grid on
% xlabel('Multipath Error (m)', 'Interpreter', 'latex', 'FontSize',16)
% ylabel('Positioning Error (m)', 'Interpreter', 'latex', 'FontSize',16)
% xlim([0, 460])
% %%
% % figure(2)
% % error_2SP_2sat = load("result/multipathErrorVSposError/pos_error_2SP_500_multipath_error_2sat.mat");
% % error_2SP_2sat = error_2SP_2sat.pos_error_2SP;
% % error_2SP_4sat = load("result/multipathErrorVSposError/pos_error_2SP_500_multipath_error.mat");
% % error_2SP_4sat = error_2SP_4sat.pos_error_2SP;
% % error_2SP_6sat = load("result/multipathErrorVSposError/pos_error_2SP_500_multipath_error_6sat.mat");
% % error_2SP_6sat = error_2SP_6sat.pos_error_2SP;
% % 
% % error_DPE_2sat = load("result/multipathErrorVSposError/pos_error_DPE_500_multipath_error_2sat.mat");
% % error_DPE_2sat = error_DPE_2sat.pos_error_DPE;
% % error_DPE_4sat = load("result/multipathErrorVSposError/pos_error_DPE_500_multipath_error.mat");
% % error_DPE_4sat = error_DPE_4sat.pos_error_DPE;
% % error_DPE_6sat = load("result/multipathErrorVSposError/pos_error_DPE_500_multipath_error_6sat.mat");
% % error_DPE_6sat = error_DPE_6sat.pos_error_DPE;
% % plot(mulipath_error_candidates, error_2SP_2sat, 'Color','#51f0f0', 'LineWidth', 2, 'LineStyle',':');
% % hold on
% % plot(mulipath_error_candidates, error_2SP_4sat, 'Color','#067bc9', 'LineWidth', 2, 'LineStyle',':');
% % plot(mulipath_error_candidates, error_2SP_6sat, 'Color', '#0000ff', 'LineWidth', 2, 'LineStyle',':');
% % 
% % plot(mulipath_error_candidates, error_DPE_2sat, 'Color','#ffb0c0', 'LineWidth', 2);
% % plot(mulipath_error_candidates, error_DPE_4sat, 'Color','#f53b79', 'LineWidth', 2);
% % plot(mulipath_error_candidates, error_DPE_6sat, 'Color', '#ba1130', 'LineWidth', 2);
% % legend("2SP N=2", "2SP N=4", "2SP N=6", "DPE N=2", "DPE N=4", "DPE N=6",...
% %     'Interpreter', 'latex', 'FontSize',14)
% % grid on
% % xlabel('Multipath Error (m)', 'Interpreter', 'latex', 'FontSize',16)
% % ylabel('Positioning Error (m)', 'Interpreter', 'latex', 'FontSize',16)
% % xlim([0, 460])
%%
% figure(2)
% error_2SP_2sat = load("result/multipathErrorVSposError/pos_error_2SP_500_multipath_error_2sat.mat");
% error_2SP_2sat = error_2SP_2sat.pos_error_2SP;
% error_2SP_4sat = load("result/multipathErrorVSposError/pos_error_2SP_500_multipath_error.mat");
% error_2SP_4sat = error_2SP_4sat.pos_error_2SP;
% 
% error_DPE_2sat = load("result/multipathErrorVSposError/pos_error_DPE_500_multipath_error_2sat.mat");
% error_DPE_2sat = error_DPE_2sat.pos_error_DPE;
% error_DPE_4sat = load("result/multipathErrorVSposError/pos_error_DPE_500_multipath_error.mat");
% error_DPE_4sat = error_DPE_4sat.pos_error_DPE;
% 
% mulipath_error_candidates = (0:10:500);
% plot(mulipath_error_candidates, error_2SP_4sat, 'Color','#0072bd', 'LineWidth', 2, 'LineStyle',':');
% hold on
% plot(mulipath_error_candidates, error_2SP_2sat, 'Color','#6bcffa', 'LineWidth', 2, 'LineStyle',':');
% 
% plot(mulipath_error_candidates, error_DPE_4sat, 'Color','#c90428', 'LineWidth', 2);
% plot(mulipath_error_candidates, error_DPE_2sat, 'Color','#ff6929', 'LineWidth', 2);
% 
% legend("2SP N=4", "2SP N=2", "DPE N=4", "DPE N=2",...
%     'Interpreter', 'latex', 'FontSize',14)
% grid on
% xlabel('Multipath Error (m)', 'Interpreter', 'latex', 'FontSize',16)
% ylabel('Positioning Error (m)', 'Interpreter', 'latex', 'FontSize',16)
% xlim([0, 460])

%% multipath error LMSCM
% figure(3)
% user_type = 'Car';
% scenario = 'Suburban';
% error_2SP_0 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_LOS.mat", user_type, scenario));
% error_2SP_0 = error_2SP_0.pos_error_2SP;
% error_2SP_1 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_3NLOS.mat", user_type, scenario));
% error_2SP_1 = error_2SP_1.pos_error_2SP;
% error_2SP_2 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_LOS_2.mat", user_type, scenario));
% error_2SP_2 = error_2SP_2.pos_error_2SP;
% 
% time = 0:0.01:4.99*3-0.01;
% error_2SP_cs = [error_2SP_0', error_2SP_1', error_2SP_2'];
% plot(error_2SP_cs)
% hold on 
% 
% 
% user_type = 'Car';
% scenario = 'Urban';
% error_2SP_0 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_LOS.mat", user_type, scenario));
% error_2SP_0 = error_2SP_0.pos_error_2SP;
% error_2SP_1 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_3NLOS.mat", user_type, scenario));
% error_2SP_1 = error_2SP_1.pos_error_2SP;
% error_2SP_2 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_LOS_2.mat", user_type, scenario));
% error_2SP_2 = error_2SP_2.pos_error_2SP;
% 
% time = 0:0.01:3*4.99-0.01;
% error_2SP_cu = [error_2SP_0', error_2SP_1', error_2SP_2'];
% plot(error_2SP_cu)
% 
% 
% user_type = 'Pedestrian';
% scenario = 'Suburban';
% error_2SP_0 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_LOS.mat", user_type, scenario));
% error_2SP_0 = error_2SP_0.pos_error_2SP;
% error_2SP_1 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_3NLOS.mat", user_type, scenario));
% error_2SP_1 = error_2SP_1.pos_error_2SP;
% error_2SP_2 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_LOS_2.mat", user_type, scenario));
% error_2SP_2 = error_2SP_2.pos_error_2SP;
% 
% time = 0:0.01:3*4.99-0.01;
% error_2SP_ps = [error_2SP_0', error_2SP_1', error_2SP_2'];
% plot(error_2SP_ps)
% 
% hold on
% user_type = 'Pedestrian';
% scenario = 'Urban';
% error_2SP_0 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_LOS.mat", user_type, scenario));
% error_2SP_0 = error_2SP_0.pos_error_2SP;
% error_2SP_1 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_3NLOS.mat", user_type, scenario));
% error_2SP_1 = error_2SP_1.pos_error_2SP;
% error_2SP_2 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_LOS_2.mat", user_type, scenario));
% error_2SP_2 = error_2SP_2.pos_error_2SP;
% time = 0:0.01:3*4.99-0.01;
% error_2SP_pu = [error_2SP_0', error_2SP_1', error_2SP_2'];
% plot(error_2SP_pu)
% legend("CS", "CU", "PS", "PU")

%%
% figure(4)
% user_type = 'Car';
% scenario = 'Suburban';
% error_DPE_0 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_LOS.mat", user_type, scenario));
% error_DPE_0 = error_DPE_0.pos_error_DPE;
% error_DPE_1 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_3NLOS.mat", user_type, scenario));
% error_DPE_1 = error_DPE_1.pos_error_DPE;
% error_DPE_2 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_LOS_2.mat", user_type, scenario));
% error_DPE_2 = error_DPE_2.pos_error_DPE;
% 
% time = 0:0.01:4.99*3-0.01;
% error_DPE_cs = [error_DPE_0', error_DPE_1', error_DPE_2'];
% plot(error_DPE_cs)
% hold on 
% 
% 
% 
% 
% user_type = 'Car';
% scenario = 'Urban';
% error_DPE_0 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_LOS.mat", user_type, scenario));
% error_DPE_0 = error_DPE_0.pos_error_DPE;
% error_DPE_1 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_3NLOS.mat", user_type, scenario));
% error_DPE_1 = error_DPE_1.pos_error_DPE;
% error_DPE_2 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_LOS_2.mat", user_type, scenario));
% error_DPE_2 = error_DPE_2.pos_error_DPE;
% 
% time = 0:0.01:4.99*3-0.01;
% error_DPE_cu = [error_DPE_0', error_DPE_1', error_DPE_2'];
% plot(error_DPE_cu)
% 
% 
% user_type = 'Pedestrian';
% scenario = 'Suburban';
% error_DPE_0 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_LOS.mat", user_type, scenario));
% error_DPE_0 = error_DPE_0.pos_error_DPE;
% error_DPE_1 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_3NLOS.mat", user_type, scenario));
% error_DPE_1 = error_DPE_1.pos_error_DPE;
% error_DPE_2 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_LOS_2.mat", user_type, scenario));
% error_DPE_2 = error_DPE_2.pos_error_DPE;
% 
% time = 0:0.01:4.99*3-0.01;
% error_DPE_ps = [error_DPE_0', error_DPE_1', error_DPE_2'];
% plot(error_DPE_ps)
% 
% 
% user_type = 'Pedestrian';
% scenario = 'Urban';
% error_DPE_0 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_LOS.mat", user_type, scenario));
% error_DPE_0 = error_DPE_0.pos_error_DPE;
% error_DPE_1 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_3NLOS.mat", user_type, scenario));
% error_DPE_1 = error_DPE_1.pos_error_DPE;
% error_DPE_2 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_LOS_2.mat", user_type, scenario));
% error_DPE_2 = error_DPE_2.pos_error_DPE;
% 
% time = 0:0.01:4.99*3-0.01;
% error_DPE_ps = [error_DPE_0', error_DPE_1', error_DPE_2'];
% plot(error_DPE_ps)
% legend("CS", "CU", "PS", "PU")

%% 
% fig_id = 5;
% for user_type = ["Car", "Pedestrian"]
%     for scenario = ["Suburban", "Urban"]
%         figure(fig_id)
%         error_2SP_0 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_LOS.mat", user_type, scenario));
%         error_2SP_0 = error_2SP_0.pos_error_2SP;
%         error_2SP_1 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_3NLOS.mat", user_type, scenario));
%         error_2SP_1 = error_2SP_1.pos_error_2SP;
%         error_2SP_2 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_LOS_2.mat", user_type, scenario));
%         error_2SP_2 = error_2SP_2.pos_error_2SP;
%         error_2SP = [error_2SP_0', error_2SP_1', error_2SP_2'];
% 
%         error_DPE_0 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_LOS.mat", user_type, scenario));
%         error_DPE_0 = error_DPE_0.pos_error_DPE;
%         error_DPE_1 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_3NLOS.mat", user_type, scenario));
%         error_DPE_1 = error_DPE_1.pos_error_DPE;
%         error_DPE_2 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_LOS_2.mat", user_type, scenario));
%         error_DPE_2 = error_DPE_2.pos_error_DPE;
%         error_DPE = [error_DPE_0', error_DPE_1', error_DPE_2'];
%         time = 0:0.01:4.99*3-0.01;
% 
%         plot(time, error_2SP, "LineWidth", 1, "Color", "#0072BD")
%         hold on
%         plot(time, error_DPE, "LineWidth", 1, "Color", "#77AC30")
%         xlabel("Navigation time (s)", 'Interpreter', 'latex', 'FontSize',16)
%         ylabel("Positioning Error (m)", 'Interpreter', 'latex', 'FontSize',16)
%         legend("2SP", "DPE", 'Interpreter', 'latex', 'FontSize',16)
%         ylim([0, 400])
%         save_name = sprintf("result/scenario/pos_error_%s_%s.fig", user_type, scenario);
%         saveas(gcf, save_name)
% 
%         fig_id = fig_id + 1;
% 
%     end
% end

%%
% fig_id = 9;
% for user_type = ["Car", "Pedestrian"]
%     for scenario = ["Suburban", "Urban"]
%         figure(fig_id)
%         error_2SP_0 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_LOS.mat", user_type, scenario));
%         error_2SP_0 = error_2SP_0.pos_error_2SP;
%         error_2SP_1 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_3NLOS.mat", user_type, scenario));
%         error_2SP_1 = error_2SP_1.pos_error_2SP;
%         error_2SP_2 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_LOS_2.mat", user_type, scenario));
%         error_2SP_2 = error_2SP_2.pos_error_2SP;
%         error_2SP = [error_2SP_0', error_2SP_1', error_2SP_2'];
% 
%         error_DPE_0 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_LOS.mat", user_type, scenario));
%         error_DPE_0 = error_DPE_0.pos_error_DPE;
%         error_DPE_1 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_3NLOS.mat", user_type, scenario));
%         error_DPE_1 = error_DPE_1.pos_error_DPE;
%         error_DPE_2 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_LOS_2.mat", user_type, scenario));
%         error_DPE_2 = error_DPE_2.pos_error_DPE;
%         error_DPE = [error_DPE_0', error_DPE_1', error_DPE_2'];
%         time = 0:0.01:4.99*3-0.01;
% 
%         plot(time, error_2SP, "LineWidth", 1, "Color", "#0072BD")
%         hold on
%         plot(time, error_DPE, "LineWidth", 1, "Color", "#77AC30")
%         xlim([1, 4.9])
%         ylim([0, 100])
%         set(gca, 'XTick', []);
%         xlabel('');
%         set(gca, 'YTick', [0, 20, 40, 60, 80, 100]);
%         ylabel('');
%         save_name = sprintf("result/scenario/pos_error_%s_%s_zoom_in.fig", user_type, scenario);
%         saveas(gcf, save_name)
% 
%         fig_id = fig_id + 1;
% 
%     end
% end

%% 
quan_80th_2SP = zeros(4, 2);
quan_50th_2SP = zeros(4, 2);
quan_80th_DPE = zeros(4, 2);
quan_50th_DPE = zeros(4, 2);
ind = 1;
for user_type = ["Car", "Pedestrian"]
    for scenario = ["Suburban", "Urban"]
        error_2SP_0 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_LOS.mat", user_type, scenario));
        error_2SP_0 = error_2SP_0.pos_error_2SP;
        error_2SP_1 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_3NLOS.mat", user_type, scenario));
        error_2SP_1 = error_2SP_1.pos_error_2SP;
        error_2SP_2 = load(sprintf("result/scenario/pos_error_2SP_%s_%s_5s_LOS_2.mat", user_type, scenario));
        error_2SP_2 = error_2SP_2.pos_error_2SP;
        error_2SP_LOS = [error_2SP_0', error_2SP_2', error_2SP_1(1:50)', error_2SP_1(451:499)'];
        error_2SP_NLOS = [error_2SP_1(51:450)'];
        quan_80th_2SP(ind, 1) = quantile(error_2SP_LOS, 0.8);
        quan_80th_2SP(ind, 2) = quantile(error_2SP_NLOS, 0.8);
        quan_50th_2SP(ind, 1) = quantile(error_2SP_LOS, 0.5);
        quan_50th_2SP(ind, 2) = quantile(error_2SP_NLOS, 0.5);
        


        error_DPE_0 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_LOS.mat", user_type, scenario));
        error_DPE_0 = error_DPE_0.pos_error_DPE;
        error_DPE_1 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_3NLOS.mat", user_type, scenario));
        error_DPE_1 = error_DPE_1.pos_error_DPE;
        error_DPE_2 = load(sprintf("result/scenario/pos_error_DPE_%s_%s_5s_LOS_2.mat", user_type, scenario));
        error_DPE_2 = error_DPE_2.pos_error_DPE;
        error_DPE_LOS = [error_DPE_0', error_DPE_2', error_DPE_1(1:50)', error_DPE_1(451:499)'];
        error_DPE_NLOS = [error_DPE_1(51:450)'];
        quan_80th_DPE(ind, 1) = quantile(error_DPE_LOS, 0.8);
        quan_80th_DPE(ind, 2) = quantile(error_DPE_NLOS, 0.8);
        quan_50th_DPE(ind, 1) = quantile(error_DPE_LOS, 0.5);
        quan_50th_DPE(ind, 2) = quantile(error_DPE_NLOS, 0.5);
        ind = ind + 1;
    end
end