function [x, y, z] = utm2ecef(E, N, U, utmZone)
    % Converts UTM coordinates to ECEF coordinates
    % Inputs:
    %   E, N - UTM Easting and Northing
    %   U - Altitude (Up coordinate)
    %   utmZone - UTM zone (e.g., '33T')
    % Outputs:
    %   x, y, z - ECEF coordinates

    % Define WGS84 constants
    a = 6378137.0; % semi-major axis
    e = 0.08181919084; % first eccentricity

    % Extract the zone number and the hemisphere
    zone = utmZone;
    hemisphere = utmZone(end);
    
    % Calculate the central meridian of the zone
    lon0 = (zone - 1) * 6 - 180 + 3;
    lon0 = deg2rad(lon0);

    % Remove the false northing for southern hemisphere
    if hemisphere == 'S'
        N = N - 10000000;
    end

    % Remove the false easting
    E = E - 500000;

    % Calculate the footpoint latitude
    k0 = 0.9996;
    e1 = (1 - sqrt(1 - e^2)) / (1 + sqrt(1 - e^2));
    M = N / k0;
    mu = M / (a * (1 - e^2 / 4 - 3 * e^4 / 64 - 5 * e^6 / 256));
    
    phi1 = mu + (3 * e1 / 2 - 27 * e1^3 / 32) * sin(2 * mu) + ...
        (21 * e1^2 / 16 - 55 * e1^4 / 32) * sin(4 * mu) + ...
        (151 * e1^3 / 96) * sin(6 * mu) + (1097 * e1^4 / 512) * sin(8 * mu);

    % Calculate latitude and longitude
    C1 = e^2 / (1 - e^2) * cos(phi1)^2;
    T1 = tan(phi1)^2;
    N1 = a / sqrt(1 - e^2 * sin(phi1)^2);
    R1 = N1 * (1 - e^2) / (1 - e^2 * sin(phi1)^2);
    D = E / (N1 * k0);

    lat = phi1 - (N1 * tan(phi1) / R1) * (D^2 / 2 - (5 + 3 * T1 + 10 * C1 - 4 * C1^2 - 9 * e^2) * D^4 / 24 + ...
        (61 + 90 * T1 + 298 * C1 + 45 * T1^2 - 252 * e^2 - 3 * C1^2) * D^6 / 720);

    lon = lon0 + (D - (1 + 2 * T1 + C1) * D^3 / 6 + ...
        (5 - 2 * C1 + 28 * T1 - 3 * C1^2 + 8 * e^2 + 24 * T1^2) * D^5 / 120) / cos(phi1);

    % Convert radians to degrees
    lat = rad2deg(lat);
    lon = rad2deg(lon);

    % Convert geodetic coordinates to ECEF coordinates
    [x, y, z] = geodetic2ecef(lat, lon, U);
end

function [x, y, z] = geodetic2ecef(lat, lon, alt)
    % Converts geodetic coordinates to ECEF coordinates
    % Inputs:
    %   lat - Latitude in degrees
    %   lon - Longitude in degrees
    %   alt - Altitude in meters
    % Outputs:
    %   x, y, z - ECEF coordinates

    % Define WGS84 constants
    a = 6378137.0; % semi-major axis
    e = 0.08181919084; % first eccentricity

    % Convert latitude and longitude to radians
    lat = deg2rad(lat);
    lon = deg2rad(lon);

    % Calculate ECEF coordinates
    N = a / sqrt(1 - e^2 * sin(lat)^2);
    x = (N + alt) * cos(lat) * cos(lon);
    y = (N + alt) * cos(lat) * sin(lon);
    z = (N * (1 - e^2) + alt) * sin(lat);
end
