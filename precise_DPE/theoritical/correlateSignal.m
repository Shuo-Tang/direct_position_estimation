function r = correlateSignal(sigen,config,x_delay_noise)


%% load configuration file
eval(config)

x_local = sigen.x_local;
NsamplesLocal = sigen.NsamplesLocal;
NsamplesData = sigen.NsamplesData;

%% memory allocation
r = cell(7,1);
if DPEflag == 0
    nfd = sigen.nfd;
    for kSV=1:numSV
        r{kSV,:} = zeros(nfd,NsamplesLocal);
    end
end
%% correlation
if DPEflag == 0
    for kSV=1:numSV
        x_delay_fft = fft(x_delay_noise,NsamplesData);
        for kfd = 1:nfd           
            x_local_fft = fft(x_local{kSV}(kfd,:),NsamplesLocal);
%             r = abs(ifft(x_delay_fft.*conj(x_local_fft))).^2;
            r{kSV,:}(kfd,:)= abs(ifft(x_delay_fft.*conj(x_local_fft))).^2;
        end
    end
    
else    
    SearchRange = sigen.SearchRange;
    for kSV=1:numSV
        for kX = 1:(2*SearchRange+1)
            for kY = 1:(2*SearchRange+1)
                x = x_local{kSV}(kX,kY,:);
                x = x(:).';
                r{kSV,:}(kX,kY)= abs(sum((x_delay_noise.*conj(x)))).^2;
            end
        end
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



