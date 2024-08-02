function PosErrDPE = DPEbruPVT(r,sigen,config)
%% load configuration file
eval(config)

SearchRange = sigen.SearchRange;
SearchDensity = sigen.SearchDensity;
UserPosition = param.UserPosition;
%% Search for max CAF
Cost = 0;
m = (0-SearchRange*SearchDensity:SearchDensity:0+SearchRange*SearchDensity);
n = (0-SearchRange*SearchDensity:SearchDensity:0+SearchRange*SearchDensity);
% figure,
for kSV=1:7
Cost = Cost + r{kSV,:}(1:(2*SearchRange+1),1:(2*SearchRange+1));
% subplot(2,4,kSV)
% mesh(m,n,r{kSV,:}(1:(2*SearchRange+1),1:(2*SearchRange+1)))
% title(strcat('CAF for Satellite ',num2str(kSV)))
% xlabel('$$y-\hat{y}$$','Interpreter','Latex')
% ylabel('$$x-\hat{x}$$','Interpreter','Latex')
% zlabel('CAF')
end
% subplot(2,4,8),mesh(m,n,Cost)
% title('Total Cost Function')
% xlabel('$$y-\hat{y}$$','Interpreter','Latex')
% ylabel('$$x-\hat{x}$$','Interpreter','Latex')
% zlabel('CAF')

[m0,n0] = find(Cost == max(Cost(:)))
ran = randperm(length(m0),1);
m0 = m0(ran);
n0 = n0(ran);
EstRxPVT = [UserPosition(1)+(SearchRange+1-m0)*SearchDensity,...
    UserPosition(2)+(SearchRange+1-n0)*SearchDensity,...
    UserPosition(3)];

PosErrDPE=norm(EstRxPVT-UserPosition);
end

