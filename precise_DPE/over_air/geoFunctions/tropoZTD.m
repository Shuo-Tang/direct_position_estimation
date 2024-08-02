function [trop_hs, trop_wet, z] = tropoZTD(pos, el, humi)
    % saastamonien tropospheric delay model
    
    temp0 = 15; % temperature at sea level
    
    if pos(3) < -100 || pos(3) > 1e4 || el <= 0
        trop_hs = 0;
        trop_wet = 0;
        z = 0;
        return;
    end
    
    hgt = max(pos(3), 0);
    
    % standard atmosphere
    pres = 1013.25 * (1 - 2.2557e-5 * hgt) ^ 5.2568;
    temp = temp0 - 6.5e-3 * hgt + 273.16;
    e = 6.108 * humi * exp((17.15 * temp - 4684.0) / (temp - 38.45));
    
    % saastamoinen model
    z = pi / 2.0 - el;
    trop_hs = 0.0022768 * pres / (1.0 - 0.00266 * cos(2 * pos(1)) - ...
        0.00028e-3 * hgt) / cos(z);
    trop_wet = 0.002277 * (1255.0 / temp + 0.05) * e / cos(z);
        
end

