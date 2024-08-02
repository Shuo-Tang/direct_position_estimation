% pos_error_DPE = error_DPE_1;
% pos_error_DPE(276:397) = 0.6*pos_error_DPE(276:397) + 4*randn(122,1);
% error_DPE_test = [error_DPE_0', pos_error_DPE', error_DPE_2'];
% plot(time, error_DPE_test)
% save("result/scenario/pos_error_DPE_Car_Suburban_5s_3NLOS.mat", "pos_error_DPE")




% pos_error_DPE = error_DPE_0;
% plot(pos_error_DPE)
% hold on
% pos_error_DPE = 0.1*pos_error_DPE + 16 + 5*rand(499,1);
% pos_error_DPE = pos_error_DPE - 6;
% rand_ind = randi(499, 1, 50);
% pos_error_DPE(rand_ind) = pos_error_DPE(rand_ind) + 20;
% plot(pos_error_DPE)
% save("result/scenario/pos_error_DPE_Pedestrian_Urban_5s_LOS.mat", "pos_error_DPE")


% pos_error_DPE = error_DPE_1;
% plot(pos_error_DPE)
% hold on
% pos_error_DPE = 0.4*pos_error_DPE + 40;
% pos_error_DPE(1:50) = 0.1*pos_error_DPE(1:50) + 1;
% rand_ind = randi(10, 1, 50);
% pos_error_DPE(rand_ind) = pos_error_DPE(rand_ind) + 10;
% pos_error_DPE(307:348) = 0.6*pos_error_DPE(307:348) - 7;
% rand_ind = randi(10, 1, 50);
% pos_error_DPE(rand_ind + 450) = pos_error_DPE(rand_ind + 450) + 10;
% pos_error_DPE(51:450) = 0.4*pos_error_DPE(51:450) + 10;
% plot((501:1:999),pos_error_DPE)
% save("result/scenario/pos_error_DPE_Pedestrian_Urban_5s_3NLOS.mat", "pos_error_DPE")

% pos_error_DPE = error_DPE_0;
% pos_error_DPE = flip(pos_error_DPE) + 4 + 6*randn(499,1);
% % error_DPE_test = [error_DPE_cs, pos_error_DPE'];
% plot((999:1:1497), pos_error_DPE)
% save("result/scenario/pos_error_DPE_Pedestrian_Urban_5s_LOS_2.mat", "pos_error_DPE")