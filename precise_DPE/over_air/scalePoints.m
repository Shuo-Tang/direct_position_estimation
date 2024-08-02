function newPoints = scalePoints(points, scaleFactor)
    % points: Nx3 matrix of 3D points (each row is a point [x, y, z])
    % scaleFactor: scaling factor (greater than 1 to spread out, between 0 and 1 to concentrate)
    
    % Calculate the centroid of the points
    centroid = mean(points);
    
    % Translate points so that centroid is at the origin
    translatedPoints = points - centroid;
    
    % Apply the scaling factor
    scaledPoints = translatedPoints * scaleFactor;
    
    % Translate points back to their original position
    newPoints = scaledPoints + centroid;
end