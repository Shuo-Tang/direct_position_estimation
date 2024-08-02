function ep = time2epoch(t)
    % convert time to epoch
    
    mday = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, ...
            30, 31, 31, 30, 31, 30, 31, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, ...
            30, 31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

    days = floor(t.time / 86400);
    sec = floor(t.time - days * 86400);
    day = mod(days, 1461);
    
    for mon = 1:48
        if day >= mday(mon)
            day = day - mday(mon);
        else
            break;
        end
    end
    
    ep = zeros(1, 6);
    ep(1) = 1970 + floor(days / 1461) * 4 + floor(mon / 12);
    ep(2) = mod(mon, 12) + 1;
    ep(3) = day + 1;
    ep(4) = floor(sec / 3600);
    ep(5) = floor(mod(sec, 3600) / 60);
    ep(6) = mod(sec, 60) + t.sec;
    end