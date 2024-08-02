function r = correlateSignal0(sigen,config,x_delay_noise)


%% load configuration file
eval(config)

x_local = sigen.x_delay;
NsamplesLocal = sigen.NsamplesLocal;
NsamplesData = sigen.NsamplesData;

%% memory allocation
r = zeros(7,1);
%% correlation
if DPEflag == 1
    for kSV=1:numSV
        r(kSV,:)= real(sum((x_delay_noise(kSV,:).*conj(x_local(kSV,:)))))^2;
    end
else
    for kSV=1:numSV
        r(kSV,:)= abs(sum((x_delay_noise(kSV,:).*conj(x_local(kSV,:)))))^2;
    end
end

%% visualization
% fdcan = (-5000:100:5000);
% delay = dt*(1:1:NsamplesLocal);
% figure, mesh(delay,fdcan,r{3,:})
% 
% x = (1:1:1001);
% y = (1:1:1001);
% figure, mesh(x,y,r{2,:})

%%
% Cost = 0;
% m = (0-SearchRange:1:0+SearchRange);
% n = (0-SearchRange:1:0+SearchRange);
% figure,
% for kSV=1:7
% Cost = Cost + r{kSV,:}(1:(2*SearchRange+1),1:(2*SearchRange+1));
% subplot(2,4,kSV)
% mesh(m,n,r{kSV,:}(1:(2*SearchRange+1),1:(2*SearchRange+1)))
% title(strcat('CAF for Satellite ',num2str(kSV)))
% xlabel('$$y-\hat{y}$$','Interpreter','Latex')
% ylabel('$$x-\hat{x}$$','Interpreter','Latex')
% zlabel('CAF')
% end
% subplot(2,4,8),mesh(m,n,Cost)
% title('Total Cost Function')
% xlabel('$$y-\hat{y}$$','Interpreter','Latex')
% ylabel('$$x-\hat{x}$$','Interpreter','Latex')
% zlabel('CAF')
% 
% [m0,n0] = find(Cost == max(Cost(:)))



