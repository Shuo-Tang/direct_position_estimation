function [x,y] = SetCandidates(UserPosition,range,density)
x = (UserPosition(1)-range*density):density:(UserPosition(1)+range*density);
y = (UserPosition(2)-range*density):density:(UserPosition(2)+range*density);
end

