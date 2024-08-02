% config = 'ConfigFile';
% eval(config)
% trueParam.UserPosition = param.UserPosition;
% trueParam.UserVelocity = param.UserVelocity;
% trueParam.deltaT = param.deltaT;
% trueParam.deltaTdot = param.deltaTdot;
% %%
% sigen0 = signalGen(config,trueParam);
% CNo = 55*ones(numSV,1);
% x_delay_noise = receivedSignal0(sigen0,config,CNo);
% r0 = correlateSignal0(sigen0,config,x_delay_noise);
% sum(r0)
% %%
% gamma1 = trueParam;
% gamma1.UserPosition(1,1) = param.UserPosition(1,1) -2.517722708173096;
% gamma1.UserPosition(1,2) = param.UserPosition(1,2) -1.463650833815336;
% sigen1 = signalGen(config,gamma1);
% r1 = correlateSignal0(sigen1,config,x_delay_noise);
% sum(r1)
% %%
% gamma2 = trueParam;
% gamma2.UserPosition(1,1) = param.UserPosition(1,1) +0.541540623176843;
% gamma2.UserPosition(1,2) = param.UserPosition(1,2) +0.527588421478868;
% sigen2 = signalGen(config,gamma2);
% r2 = correlateSignal0(sigen2,config,x_delay_noise);
% sum(r2)
% %%
% dX = (-5:0.1:5);
% dY = (-5:0.1:5);
% for nx = 1:length(dX)
%     for ny = 1:length(dY)
%         gamma = trueParam;
%         gamma.UserPosition(1,1) = param.UserPosition(1,1) +dX(nx);
%         gamma.UserPosition(1,2) = param.UserPosition(1,2) +dY(ny);
%         sigen = signalGen(config,gamma);
%         r_grid(nx,ny) = sum(correlateSignal0(sigen,config,x_delay_noise));
%         fprintf("grid coomputation %d, %d\n",nx,ny)
%     end
% end
% figure,
% [X,Y] = meshgrid(-5:0.1:5);
% mesh(X,Y,r_grid)
%%
M = 5;
P = 1e-7*randn(5,3);
xi = randn*5*ones(M,1);
Xi = diag(xi);
gamma = randn*10*ones(M,1);
Gamma = diag(gamma);
f = 1.5e9;

A = 2*P'*Xi*Gamma*P + 8*pi^2*f^2*P'*Gamma*P;
B = 2*P'*Xi*gamma + 8*pi^2*f^2*P'*gamma;
C = B';
D = 2*xi'*gamma + 8*pi*f^2*sum(gamma);
I = [A B;C D]
bound = inv(I)
inv(A)
