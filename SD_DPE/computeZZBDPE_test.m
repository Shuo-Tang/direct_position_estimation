function fZZLB_DPE = computeZZBDPE_test(config,sigen)
%% load data
eval(config)

NsamplesData = sigen.NsamplesData;
x_delay = sigen.x_delay;
s = sigen.s;
UserPosition = param.UserPosition;
SatPosition = param.SatPosition;
deltaTdot = param.deltaTdot;
Tchip = CodePeriod / CodeLen;
dt = 1/fs;
B_2=sum((diff(s(1,:))/dt).^2)/sum(s(1,:).^2);

%% prepare some intermediate variables
fZZLB_DPE = zeros(1,length(CNosim));
fZZLB_CB = zeros(1,length(CNosim));
%% time invariant
normPos = vecnorm(param.UserPosition-param.SatPosition,2,2); % 3d norm
P= 1/c*(UserPosition-SatPosition)./normPos;

%% comput ZZB
T=CodePeriod;
M=numSV;
SNRdb=CNosim0+10*log10(T);
D=T*c;
Rd=D^2/12*eye(2);
priorCov = (RangeSet.*DensitySet).^2;


for SNRloop=1:length(SNRdb)
    SNR=10^(SNRdb(SNRloop)/10);
%     alpha = sqrt(10^(SNRdb(SNRloop)/10)/NsamplesData);
    I11 = 2*P'*B_2*SNR*P + 8*pi^2*f0^2*P'*SNR*P;
    I22 = 2*M*B_2*SNR + 8*pi^2*f0^2*M*SNR;
    I12 = 2*P'*B_2*SNR*ones(M,1) + 8*pi^2*f0^2*P'*SNR*ones(M,1);
    I21 = I12';
    IFIM = [I11,I12;I21,I22];
    
    ZZLB_DPE = priorCov(SNRloop)*2*q(sqrt(M*SNR))+(IFIM)^(-1)*gammainc(M*SNR/2,3/2);
    fZZLB_DPE(SNRloop)= sqrt(ZZLB_DPE(1,1)+ZZLB_DPE(2,2)+ZZLB_DPE(3,3));
    fZZLB_CB(SNRloop) = sqrt(ZZLB_DPE(4,4));
end
end

