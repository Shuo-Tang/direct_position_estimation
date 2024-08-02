function [localMax, cost, path] = generalizedPatternSearch2d(usrInfo, satInfo, rawSignal, caCodesTable,...
        timeDiff, searchD, stepSize, stepDecay, convergence)
%GENERALIZEDPATTERNSEARCH This function is for generalized pattern search
%optimization method in 2D scenario
%   Detailed explanation goes here
step0 = stepSize;
r = generateCAFs(usrInfo, satInfo, rawSignal, caCodesTable, timeDiff);
path = [usrInfo(1:2), r];
while stepSize > convergence
    improved  = false;
    for i = 1:size(searchD,2)
        usrInfo1 = usrInfo;
        usrInfo1(1:2) = usrInfo1(1:2) + stepSize*searchD(:,i)';
        r1 = generateCAFs(usrInfo1, satInfo, rawSignal, caCodesTable, timeDiff);
        if r1 > r
            usrInfo = usrInfo1;
            r = r1;
            path = [path; usrInfo(1:2), r];
            improved = true;
            stepSize = step0;
            D = searchD(:,i);
            searchD(:,i) = [];
            searchD = [D, searchD];
            break
        end
    end
    if ~improved
        stepSize = stepDecay*stepSize;
    end
end
localMax = usrInfo;
cost = r;
end


