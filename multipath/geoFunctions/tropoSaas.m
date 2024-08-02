function ddr = tropoSaas(pos, el)
%TROPO  Calculation of tropospheric correction.
%       The range correction ddr in m is to be subtracted from
%       pseudo-ranges and carrier phases
%
%ddr = tropo(sinel, hsta, p, tkel, hum, hp, htkel, hhum);
%
%   Inputs:
%       sinel   - sin of elevation angle of satellite
%       hsta    - height of station in km
%       p       - atmospheric pressure in mb at height hp
%       tkel    - surface temperature in degrees Kelvin at height htkel
%       hum     - humidity in % at height hhum
%       hp      - height of pressure measurement in km
%       htkel   - height of temperature measurement in km
%       hhum    - height of humidity measurement in km
%
%   Outputs:
%       ddr     - range correction (meters)
%
    humi = 0.7;
    temp0 = 15; % temparature at sea level

    posLLA = ecef2lla(pos);
    hgt = egm96geoid(posLLA(1), posLLA(2));

    % standard atmosphere
    pres = 1013.25 * (1 - 2.2557e-5 * hgt) ^ 5.2568;
    temp = temp0 - 6.5e-3 * hgt + 273.16;
    e = 6.108 * humi * exp((17.15 * temp - 4684.0) / (temp - 38.45));
    % saastamoinen model
    z = pi / 2 - deg2rad(el);
    ddr = 0.002277 / cos(z) * (pres + (1255/temp + 0.05) * e - tan(z)^2);
end
%%%%%%%%% end tropo.m  %%%%%%%%%%%%%%%%%%%

