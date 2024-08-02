function plot2DCAF(xCan, yCan, r, truePos)
%PLOT2DCAF plot 2d CAF value and compute the local maximum
%   
figure(100)
[x, y] = meshgrid(xCan-truePos(1), yCan-truePos(2));
% mesh(x, y, r')
h = surf(x, y, r');
h.EdgeColor = [0.7 0.7 0.7];
alpha(h, 0.5);
hold on
xlabel("Y - Y_{true} (ECEF)"); ylabel("X - X_{true} (ECEF)"); zlabel("Cost")
% plot max
rMax = max(r(:));
[yMax,xMax] = find(r == rMax);
numMax = length(yMax);
% plot3(y(yMax,1), x(1, xMax), rMax*1.1*ones(1, numMax))
% pos_max = [xCan(xMax);yCan(yMax)]';
% pos = mean(pos_max, 1);
% h_max = plot3(pos(1)-truePos(1), pos(2)-truePos(2), 1.01*rMax, "v", 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% legend(h_max, 'DPE solution');
x_max = [-2,-1,0,0,1,1,1,1,1,2,2,2,2,3,3,2,2,1,1,1,0,0,-1,-1,-2,-2,-3,-3,-3,-3,-3,-2];
y_max = [3,3,3,2,2,1,0,-1,-1,-2,-3,-4,-4,-5,-6,-6,-7,-7,-6,-5,-5,-4,-4,-3,-3,-2,-2,-1,-0,-1,2,2];
h_max = plot3(x_max, y_max, rMax*1.01*ones(1, 32), 'Color',"#A2142F","LineWidth",3, 'LineStyle',':' )
legend(h_max, 'Maximum Cost Area')
xlim([-25, 25])
ylim([-25, 25])
end

