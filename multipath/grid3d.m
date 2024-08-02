function [x, y, z] = grid3d(cen, range, density)
% generate 2d search grid

x0 = cen(1);
y0 = cen(2);
z0 = cen(3);

x = (x0 - range: density: x0 + range - 1);
y = (y0 - range: density: y0 + range - 1);
z = (z0 - 0: density: z0 + 0);
end

