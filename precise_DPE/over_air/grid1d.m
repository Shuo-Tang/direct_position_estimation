function [x] = grid1d(cen, range, density)
% generate 2d search grid

x0 = cen(1);

x = (x0 - range: density: x0 + range);
end

