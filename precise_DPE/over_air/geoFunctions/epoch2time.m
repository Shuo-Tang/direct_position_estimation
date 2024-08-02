function time = epoch2time(ep)
    % calculate time from epoch
    doy = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335];
    time = gtime_t();
    year = floor(ep(1));
    mon = floor(ep(2));
    day = floor(ep(3));

    if year < 1970 || year > 2099 || mon < 1 || mon > 12
        return;
    end
    
    days = (year - 1970) * 365 + floor((year - 1969) / 4) + doy(mon) + day - 2;
    
    if mod(year, 4) == 0 && mon >= 3
        days = days + 1;
    end
    
    sec = floor(ep(6));
    time.time = days * 86400 + floor(ep(4)) * 3600 + floor(ep(5)) * 60 + sec;
    time.sec = ep(6) - sec;
    end