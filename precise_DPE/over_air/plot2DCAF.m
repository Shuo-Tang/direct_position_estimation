function plot2DCAF(xCan, yCan, r, truePos)
%PLOT2DCAF plot 2d CAF value and compute the local maximum
%   
figure(100)
[x, y] = meshgrid(xCan-truePos(1), yCan-truePos(2));
mesh(x, y, r')
xlabel("Y - Y_{true} (ECEF)"); ylabel("X - X_{true} (ECEF)"); zlabel("Cost")
saveas(gcf,'result/real_caf_20ms_large')
[trueLLA(1), trueLLA(2), trueLLA(3)] = cart2geo(truePos(1), truePos(2), truePos(3), 5);
utmZone = findUtmZone(trueLLA(1), trueLLA(2));

% [trueE, trueN, trueU] = cart2utm(truePos(1), truePos(2), truePos(3), utmZone);
% 
% 
% [x, y] = meshgrid(xCan, yCan);
% E = zeros(size(x));
% N = zeros(size(x));
% for i = 1:length(xCan)
%     for j = 1:length(yCan)
%         [E(i, j), N(i, j), ~] = cart2utm(x(i, j), y(i, j), truePos(3), utmZone);
%         E(i, j) = E(i, j) - trueE;
%         N(i, j) = N(i, j) - trueN;
%     end
% end
% figure(101)
% mesh(E, N, r')
% xlabel("N - N_{true} (ECEF)"); ylabel("E - E_{true} (ECEF)"); zlabel("Cost")

end

