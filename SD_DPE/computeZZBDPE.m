function fZZLB_DPE = computeZZBDPE(config,sigen)
%% load data
eval(config)

NsamplesData = sigen.NsamplesData;
x_delay = sigen.x_delay;
s = sigen.s;
% s0 = sigen.s0;
UserPosition = param.UserPosition(:,1:2);
SatPosition = param.SatPosition(:,1:2);
UserVelocity = param.UserVelocity(:,1:2);
SatVelocity = param.SatVelocity(:,1:2);
deltaTdot = param.deltaTdot;
Tchip = CodePeriod / CodeLen;
dt = 1/fs;

%% prepare some intermediate variables
%% time invariant
normPos = vecnorm(param.UserPosition-param.SatPosition,2,2); % 3d norm
DiffV = repmat(UserVelocity,7,1)-SatVelocity;
dV = repmat(param.UserVelocity,7,1)-param.SatVelocity;
% DiffVx = DiffV(:,1);
% DiffVy = DiffV(:,2);
dP = repmat(param.UserPosition,7,1) - param.SatPosition;
DiffP = repmat(UserPosition,7,1) - SatPosition;
DiffPx = DiffP(:,1);
DiffPy = DiffP(:,2);
P= 1/c*(UserPosition-SatPosition)./normPos;
Q = zeros(numSV,2);
U = zeros(2,2,numSV); % U is P_dot
% W = zeros(2,2,numSV); % W is Q_dot
for kSV = 1:numSV
    normp = normPos(kSV,:);
    dp = DiffP(kSV,:);
%     dv = DiffV(kSV,:);
    dx = DiffPx(kSV,:);
    dy = DiffPy(kSV,:);
%     Vx = DiffVx(kSV,:);
%     Vy = DiffVy(kSV,:);
    
    Q(kSV,:) = -dV(kSV,:)*([eye(2);0 0]/normp-dP(kSV,:)'*dp/normp^3)*f0/c*(1+deltaTdot);
    U(1,1,kSV) = (1/normp-dx^2/normp^3)/c;
    U(1,2,kSV) = (-dy*dx)/normp^3/c;
    U(2,1,kSV) = (-dy*dx)/normp^3/c;
    U(2,2,kSV) = (1/normp-dy^2/normp^3)/c;
%     W(1,1,kSV) = Vx*(3*dx*normp^3-3*dx^3*normp)+Vy*(dy*normp^3-3*dx^2*dy*normp);
%     W(1,2,kSV) = Vx*(dy*normp^3-3*dx^2*dy*normp)+Vy*(dx*normp^3-3*dx*dy^2*normp);
%     W(2,1,kSV) = Vx*(dy*normp^3-3*dx^2*dy*normp)+Vy*(dx*normp^3-3*dx*dy^2*normp);
%     W(2,2,kSV) = Vx*(dx*normp^3-3*dx*dy^2*normp)+Vy*(3*dy*normp^3-3*dy^3*normp);
%     W(:,:,kSV) = W(:,:,kSV)./normp^6;
end
%% time variant
s_dot = diff(s,1,2)./dt;
s_dot = [zeros(numSV,1),s_dot];

% Beta = zeros(numSV,2,NsamplesData);
ifim0 = zeros(2,2,NsamplesData);
ifim1 = zeros(2,2,NsamplesData);
ifim2 = zeros(2,2,NsamplesData);
% ifim1 = zeros(2,2,NsamplesData);
fZZLB_DPE = zeros(1,length(CNosim));
% fZZLB_DPE1 = zeros(1,length(CNosim));

% for n = 1:NsamplesData
%     for kSV = 1:numSV
%         Beta(kSV,:,n) = -s_dot(kSV,n)/s(kSV,n)*P(kSV,:)+...
%             1i*2*pi*n*dt*Q(kSV,:)+1i*2*pi*f0*P(kSV,:);
% %         Beta(kSV,:,n) = -s_dot(kSV,n)/s(kSV,n)*P(kSV,:);
%     end
% end
%% comput ZZB
T=CodePeriod;
M=numSV;
SNRdb=CNosim0+10*log10(T);
D=T*c;
Rd=D^2/12*eye(2);
priorCov = (RangeSet.*DensitySet).^2;


for SNRloop=1:length(SNRdb)
    SNR=10^(SNRdb(SNRloop)/10);
    alpha = sqrt(10^(SNRdb(SNRloop)/10)/NsamplesData);
    for n = 1:NsamplesData
        H0 = zeros(2,2,numSV);
        H1 = zeros(2,2,numSV);
        H2 = zeros(2,2,numSV);
%         H1 = zeros(2,2,numSV);
        t = n*dt;
        for kSV = 1:numSV
            H2(:,:,kSV) = alpha^2*(P(kSV,:)'*s_dot(kSV,n)^2*P(kSV,:)+s(kSV,n)*s_dot(kSV,n)*U(:,:,kSV))+...
                s(kSV,n)^2*alpha^2*(Q(kSV,:)'*2*pi*t+P(kSV,:)'*2*pi*f0)*...
                (2*pi*t*Q(kSV,:)+2*pi*f0*P(kSV,:));% actually s(kSV,n)*s_dot(kSV,n)*U(:,:,kSV) doesn't matter
            H1(:,:,kSV) = alpha^2*(P(kSV,:)'*s_dot(kSV,n)^2*P(kSV,:)+s(kSV,n)*s_dot(kSV,n)*U(:,:,kSV))+...
                s(kSV,n)^2*alpha^2*(Q(kSV,:)'*2*pi*t)*(2*pi*t*Q(kSV,:));
            H0(:,:,kSV) = alpha^2*(P(kSV,:)'*s_dot(kSV,n)^2*P(kSV,:)+s(kSV,n)*s_dot(kSV,n)*U(:,:,kSV));
        end
        ifim0(:,:,n) = 2*sum(H0,3);
        ifim1(:,:,n) = 2*sum(H1,3);
        ifim2(:,:,n) = 2*sum(H2,3);
    end
    switch DPEcase
        case 0
            IFIM = sum(ifim0,3);
        case 1
            IFIM = sum(ifim1,3);
        case 2
            IFIM = sum(ifim2,3);
    end
    ZZLB_DPE= priorCov(SNRloop)*2*q(sqrt(M*SNR))+inv(IFIM)*gammainc(M*SNR/2,3/2);
    fZZLB_DPE(SNRloop)=sqrt(ZZLB_DPE(1,1)+ZZLB_DPE(2,2));
end

%% test with previous version
% B_2=sum((diff(s(1,:))/dt).^2)/sum(s(1,:).^2);
% B_2 = zeros(7,1);
% for kSV = 1:numSV
%     B_2(kSV,1) = sum((diff(s(kSV,:))/dt).^2)/sum(s(kSV,:).^2);
% end
% for SNRloop=1:length(SNRdb)
%     SNR=10^(SNRdb(SNRloop)/10);
% %     Jtau=diag(repmat(SNR*B_2*2,1,M));
%     Jtau=diag(SNR.*B_2*2);
%     J= P'*Jtau*P;
% 
%     ZZLB_DPE0= 1/16*Rd*2*q(sqrt(M*SNR))+inv(J)*gammainc(M*SNR/2,3/2);
% end

%% test with brute force
% fZZLB_DPE = zeros(1,length(SNRdb));
% 
% fn=2e6;
% order=36;
% wn=pi*fn/(pi*fs/2);
% hfilter=fir1(order,wn);
% for SNRloop=1:length(SNRdb)
% 
%     A = sqrt(10^(SNRdb(SNRloop)/10)/NsamplesData)*eye(M);
%     hcan = (1:1:100);
%     maxPOE = zeros(1,length(hcan));
%     for nh = 1:length(hcan)
%         h= hcan(nh);
%         numTest = 20;
%         deltaCan = h*[(0:1:numTest);(numTest:-1:0)]'./numTest;
%         sTheta0 = x_delay;
%         POE = zeros(1,numTest);
%         for ndelta = 1:numTest
%             sTheta0Delta = zeros(numSV,NsamplesData);
%             UserPositionC = [UserPosition+deltaCan(ndelta,:),param.UserPosition(3)];
%             for kSV = 1:numSV
%                 Code = genCAcode(SatPRN(kSV));
%                 Tchip = CodePeriod / length(Code);
%                 GeoRangeC = norm(param.SatPosition(kSV,:) - UserPositionC);
%                 RangeC = GeoRangeC + c*(param.deltaT - param.deltat(kSV));
%                 FracDelayC = mod(RangeC/c,CodePeriod);
%                 PrevNCOIndexC = -FracDelayC/Tc;
%                 ii = 1:NsamplesData;
%                 x_shiftedC = Code((1 + mod(round(PrevNCOIndexC+(ii - 1)/fs/Tchip), length(Code))));
%                 x0 = x_shiftedC;
%                 x0 = filtfilt(hfilter,1,x0);
%                 sTheta0Delta(kSV,:) = x0*sqrt((NsamplesData/sum(x0.^2)));
%             end
%             Sn = zeros(1,NsamplesData);
%             for n = 1:NsamplesData
%                 Sn(:,n) = sTheta0(:,n)'*A^2*(sTheta0(:,n)-sTheta0Delta(:,n));
%             end
%             S = sum(Sn);
%             POE(1,ndelta) = q(sqrt(S/2));
%         end
%         maxPOE(1,nh) = max(POE);
%     end
%     cov = sum(hcan.*maxPOE);
%     fZZLB_DPE(1,SNRloop) = sqrt(cov);  
% end
%% test with brute force a = [1 0]
% fZZLB_DPE = zeros(1,length(SNRdb));
% 
% fn=2e6;
% order=36;
% wn=pi*fn/(pi*fs/2);
% hfilter=fir1(order,wn);
% for SNRloop=1:length(SNRdb)
% 
%     A = sqrt(10^(SNRdb(SNRloop)/10)/NsamplesData)*eye(M);
%     hcan = (1:1:100);
%     %% for x
%     maxPOE = zeros(1,length(hcan));
%     for nh = 1:length(hcan)
%         h= hcan(nh);
%         sTheta0 = x_delay;
%         sTheta0Delta = zeros(numSV,NsamplesData);
% 
%         delta = [1 0].*h;
%         UserPositionC = [UserPosition+delta,param.UserPosition(3)];      
%         for kSV = 1:numSV
%             Code = genCAcode(SatPRN(kSV));
%             Tchip = CodePeriod / length(Code);
%             GeoRangeC = norm(param.SatPosition(kSV,:) - UserPositionC);
%             RangeC = GeoRangeC + c*(param.deltaT - param.deltat(kSV));
%             FracDelayC = mod(RangeC/c,CodePeriod);
%             PrevNCOIndexC = -FracDelayC/Tc;
%             ii = 1:NsamplesData;
%             x_shiftedC = Code((1 + mod(round(PrevNCOIndexC+(ii - 1)/fs/Tchip), length(Code))));
%             x0 = x_shiftedC;
%             x0 = filtfilt(hfilter,1,x0);
%             sTheta0Delta(kSV,:) = x0*sqrt((NsamplesData/sum(x0.^2)));
%         end
%         Sn = zeros(1,NsamplesData);
%         for n = 1:NsamplesData
%             Sn(:,n) = sTheta0(:,n)'*A^2*(sTheta0(:,n)-sTheta0Delta(:,n));
%         end
%         S = sum(Sn);
%         POE = q(sqrt(S/2));
%         maxPOE(1,nh) = POE;
%     end
%     xcov = sum(hcan.*maxPOE);
%     xZZLB_DPE(1,SNRloop) = sqrt(xcov);  
%     %% for y
%     maxPOE = zeros(1,length(hcan));
%     for nh = 1:length(hcan)
%         h= hcan(nh);
%         sTheta0 = x_delay;
%         sTheta0Delta = zeros(numSV,NsamplesData);
% 
%         delta = [0 1].*h;
%         UserPositionC = [UserPosition+delta,param.UserPosition(3)];      
%         for kSV = 1:numSV
%             Code = genCAcode(SatPRN(kSV));
%             Tchip = CodePeriod / length(Code);
%             GeoRangeC = norm(param.SatPosition(kSV,:) - UserPositionC);
%             RangeC = GeoRangeC + c*(param.deltaT - param.deltat(kSV));
%             FracDelayC = mod(RangeC/c,CodePeriod);
%             PrevNCOIndexC = -FracDelayC/Tc;
%             ii = 1:NsamplesData;
%             x_shiftedC = Code((1 + mod(round(PrevNCOIndexC+(ii - 1)/fs/Tchip), length(Code))));
%             x0 = x_shiftedC;
%             x0 = filtfilt(hfilter,1,x0);
%             sTheta0Delta(kSV,:) = x0*sqrt((NsamplesData/sum(x0.^2)));
%         end
%         Sn = zeros(1,NsamplesData);
%         for n = 1:NsamplesData
%             Sn(:,n) = sTheta0(:,n)'*A^2*(sTheta0(:,n)-sTheta0Delta(:,n));
%         end
%         S = sum(Sn);
%         POE = q(sqrt(S/2));
%         maxPOE(1,nh) = POE;
%     end
%     ycov = sum(hcan.*maxPOE);
%     yZZLB_DPE(1,SNRloop) = sqrt(ycov);
% 
%     fZZLB_DPE(1,SNRloop) = sqrt(xcov+ycov);
% end
end

