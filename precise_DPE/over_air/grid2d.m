function [x, y] = grid2d(cen, range, density)
% generate 2d search grid

x0 = cen(1);
y0 = cen(2);

x = (x0 - range: density: x0 + range);
y = (y0 - range: density: y0 + range);
end

