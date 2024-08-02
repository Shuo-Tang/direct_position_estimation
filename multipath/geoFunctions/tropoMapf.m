function [mapfh, mapfw] = tropoMapf(t, pos, el)
    % tropospheric mapping function Neil (NMF)
    nmf_coef = [
    1.2769934e-3, 1.2683230e-3, 1.2465397e-3, 1.2196049e-3, 1.2045996e-3;
    2.9153695e-3, 2.9152299e-3, 2.9288445e-3, 2.9022565e-3, 2.9024912e-3;
    62.610505e-3, 62.837393e-3, 63.721774e-3, 63.824265e-3, 64.258455e-3;
    0.0, 1.2709626e-5, 2.6523662e-5, 3.4000452e-5, 4.1202191e-5;
    0.0, 2.1414979e-5, 3.0160779e-5, 7.2562722e-5, 11.723375e-5;
    0.0, 9.0128400e-5, 4.3497037e-5, 84.795348e-5, 170.37206e-5;
    5.8021897e-4, 5.6794847e-4, 5.8118019e-4, 5.9727542e-4, 6.1641693e-4;
    1.4275268e-3, 1.5138625e-3, 1.4572752e-3, 1.5007428e-3, 1.7599082e-3;
    4.3472961e-2, 4.6729510e-2, 4.3908931e-2, 4.4626982e-2, 5.4736038e-2];
    nmf_aht = [2.53E-5, 5.49E-3, 1.14E-3]; % height correction

    if pos(3) < -1e3 || pos(3) > 20e3 || el <= 0.0
        mapfh = 0.0;
        mapfw = 0.0;
        return;
    end

    aht = nmf_aht;
    lat = rad2deg(pos(1));
    
    % year from doy 28, add half a year for southern latitudes
    y = (time2doy(t) - 28.0) / 365.25;
    y = y + 0.5 * (lat < 0);
    cosy = cos(2.0 * pi * y);
    c = interpc(nmf_coef, abs(lat));
    ah = c(1:3) - c(4:6) * cosy;
    aw = c(7:9);
    
    % ellipsoidal height is used instead of height above sea level
    dm = (1.0 / sin(el) - mapf(el, aht(1), aht(2), aht(3))) * pos(3) * 1e-3;
    mapfh = mapf(el, ah(1), ah(2), ah(3)) + dm;
    mapfw = mapf(el, aw(1), aw(2), aw(3));

    function result = mapf(el, a, b, c)
    % simple tropospheric mapping function
    sinel = sin(el);
    result = (1.0 + a / (1.0 + b / (1.0 + c))) / (sinel + (a / (sinel + b / (sinel + c))));
end
    
    

    
end
