function [x, y, z, dt] = grid4d(cen, range, density)
% generate 2d search grid

x0 = cen(1);
y0 = cen(2);
z0 = cen(3);
dt0 = cen(4);

rangeX = range(1);
rangeY = range(2);
rangeZ = range(3);
rangeDT = range(4);

densityX = density(1);
densityY = density(2);
densityZ = density(3);
densityDT = density(4);

x = (x0 - rangeX: densityX: x0 + rangeX);
y = (y0 - rangeY: densityY: y0 + rangeY);
z = (z0 - rangeZ: densityZ: z0 + rangeZ);
dt = (dt0 - rangeDT: densityDT: dt0 + rangeDT);

end

