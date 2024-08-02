clear all;
close all;
clc;

dbstop if error

plot_flag    =  0;
errors_flag  =  0;

GPS_flag     =  1;
GLO_flag     =  0;


%% 

load('TAQUIOMETER_Reference.mat')

[constellations] = goGNSS.initConstellation(GPS_flag,GLO_flag,0,0,0,0);
nSatTot = constellations.nEnabledSat;
n_sys = 1; % only GPS 


if GPS_flag
    [C1C, L1C, ~, ~, DOP, ~, SNR, ~, time_GPS, time_R, week_R, date_R, pos_R, interval, ~, ~, ~] = load_RINEX_obs('mitt084o.14o', constellations);
    [Eph, iono] = load_RINEX_nav('brdc0840.14n', constellations, 0);
end

if GLO_flag
    [C1C, L1C, ~, ~, DOP, ~, SNR, ~, time_GPS, time_R, week_R, date_R, pos_R, interval, ~, ~, ~] = load_RINEX_obs('mitt084o.14o', constellations);
    [Eph, iono] = load_RINEX_nav('brdc0840.14g', constellations, 0);
end

% SP3 = load_SP3('igs', time_GPS, week_R, constellations);

nEpochs = length(time_R);        

[phiR, lamR, hR] = cart2geod(pos_R(1), pos_R(2), pos_R(3));
phiR = phiR * 180 / pi;
lamR = lamR * 180 / pi;


%% 

dtR = zeros(nEpochs,1);      % receiver clock error
dtRdot = zeros(nEpochs-1,1); % receiver clock drift    

min_nsat_LS = 3 + n_sys;

for i = 1 : 1
    
    sat0 = find(C1C(:,i) ~= 0);

    Eph_t = rt_find_eph (Eph, time_GPS(i), nSatTot);
    [XS, dtS, XS_tx, VS_tx, time_tx, no_eph, sys] = satellite_positions(time_GPS(i), C1C(sat0,i), sat0, Eph_t, [], [], zeros(nSatTot,1), zeros(nSatTot,1), 0);
       
    [az, el, dist] = topocent(pos_R, XS);
    err_tropo = tropo_error_correction(el,hR*ones(length(el),1));
    err_iono  = iono_error_correction(phiR, lamR, az, el, time_GPS(i), iono, []);
     
    index = find(no_eph == 0);
    
    sat  = sat0(index);
    c1cA = C1C(sat0,i);
    c1c  = c1cA(index);
    % snr  = snr(index);
    el   = el(index);
    az   = az(index);
    dist = dist(index);
    XS   = XS(index,:);
    dtS  = dtS(index);
    sys  = sys(index);
    nsat = size(c1c,1);
    err_tropo = err_tropo(index);
    err_iono = err_iono(index);
    
    n = length(index);
    
    if (length(sat0) >= min_nsat_LS)
 
        for ii=0:1:5
            XR_mat = pos_R(:,ones(n,1))';
            distR_approx = sqrt(sum((XS-XR_mat).^2 ,2));
    
            A = [(pos_R(1) - XS(:,1)) ./ distR_approx, (pos_R(2) - XS(:,2)) ./ distR_approx, (pos_R(3) - XS(:,3)) ./ distR_approx, ones(n,1)];
            b = distR_approx - goGNSS.V_LIGHT*dtS + err_tropo + err_iono;
            y0 = c1c;
            
            N = (A'*A);
            if cond(N) < 100
                x = (N^-1)*A'*(y0-b);
                y_hat = A*x + b;
                v_hat = y0 - y_hat;
                pos_R  = pos_R + x(1:3);
                dtR = x(4) / goGNSS.V_LIGHT;
            end
        end
        
        XR(i,1:3) = pos_R';
    end
    
end

%%
% for j = 1 : length(XR)
%     [PHI(j), LAM(j), H(j)] = cart2geod(XR(j,1), XR(j,2), XR(j,3));
%     [EAST(j), NORTH(j), ~, ~] = cart2plan(XR(j,1), XR(j,2), XR(j,3));
% end
% 
% PHI = PHI * 180 / pi;
% LAM = LAM * 180 / pi;

% figure(1), plot(detrend(EAST),detrend(NORTH),'+')

% figure(2), plot(sqrt(detrend(NORTH).^2 + detrend(NORTH).^2))


%% PLOT

% if plot_flag
%     figure;
%     plot(LAM(3585:end),PHI(3585:end))
%     hold on;
%     plot(TAQ_Position_ELL(:,2),TAQ_Position_ELL(:,1))
%     plot_google_map 
% end



%% ERRORS

% time_sync = 3601;
% 
% EAST_syn    =   EAST (time_sync:end)';
% NORTH_syn   =   NORTH (time_sync:end)';
% 
% if errors_flag
%     for i=1:(size(TAQ_Position_ELL,1)-17)
%         error_East(i)     =     abs(TAQ_Position_UTM(i,1)-EAST_syn(i));
%         error_North(i)    =     abs(TAQ_Position_UTM(i,2)-NORTH_syn(i));
%         error(i)          =     sqrt((error_East(i)^2)+(error_North(i)^2));
%     end
% end
% 
% disp('----------------------------------')
% 
% disp('EAST: mean error')
% disp(mean(error_East))
% 
% disp('NORTH: mean error')
% disp(mean(error_North))
% 
% disp('Error: mean')
% disp (mean(error))

    