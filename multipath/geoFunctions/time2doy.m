function doy = time2doy(t)
    % convert time to day of year (DOY)
    
    ep = time2epoch(t);
    ep(2:3) = 1.0;
    ep(4:6) = 0.0;
    
    doy = (timediff(t, epoch2time(ep)) / 86400) + 1;
    
    function dt = timediff(t1, t2)
    % return time difference
    dt = t1.time - t2.time + (t1.sec - t2.sec);
end
end

    

