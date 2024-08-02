function result = interpc(coef, lat)
    % linear interpolation (lat step=15)
    
    i = floor(lat / 15.0);
    
    if i < 1
        result = coef(:, 1);
    elseif i > 4
        result = coef(:, 5);
    else
        d = lat / 15.0 - i;
        result = coef(:, i) * (1.0 - d) + coef(:, i+1) * d;
    end
end


