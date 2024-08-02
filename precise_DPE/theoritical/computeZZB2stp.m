function fZZLB_2SP = computeZZB2stp(config,sigen)
%% load data
eval(config)

nfd = sigen.nfd;
x_local = sigen.x_local{1}((nfd+1)/2,:);
UserPosition = param.UserPosition;
SatPosition = param.SatPosition;
%% compute ZZB
B_2=sum((diff(x_local(1,:))/dt).^2)/sum(x_local(1,:).^2);
T=CodePeriod;
M=numSV;

P= 1/c*(UserPosition-SatPosition)./sqrt(sum((UserPosition-SatPosition).^2,2));
varT=T^2/12;

SNRdb=CNosim+10*log10(T);
fZZLB_2SP = zeros(1,length(CNosim));

for SNRloop=1:length(SNRdb) 
    SNR=10^(SNRdb(SNRloop)/10);
    ZZLB_tau=eye(M)* (1/8*varT*2*q(sqrt(SNR))+(1/(SNR*B_2*2))*gammainc(SNR/2,3/2));
    ZZLB_2SP=(P'*P)^-1*P'*(ZZLB_tau)*((P'*P)^-1*P')';    
    fZZLB_2SP(SNRloop)=sqrt(ZZLB_2SP(1,1)+ZZLB_2SP(2,2)+ZZLB_2SP(3,3));
end
end

