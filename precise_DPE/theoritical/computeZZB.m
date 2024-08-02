x_local = sigen.x_local;

B_2=sum((diff(x_local(1,:))/dt).^2)/sum(x_local(1,:).^2);
T=CodePeriod;
D=T*c;
M=numSV;

P= 1/c*(UserPosition-SatPosition)./sqrt(sum((UserPosition-SatPosition).^2,2));

Rd=[D^2/12 0 0; 0 D^2/12 0; 0 0 D^2/12 ];
RT=[T^2/12 0 0; 0 T^2/12 0; 0 0 T^2/12 ];
varT=T^2/12;

SNRdb=CNosim+10*log10(T);
fZZLB_DPE = zeros(1,length(CNosim));
fZZLB_2SP = zeros(1,length(CNosim));

for SNRloop=1:length(SNRdb)
    
    SNR=10^(SNRdb(SNRloop)/10);
    Jtau=diag(repmat(SNR*B_2*2,1,M));
    J= P'*Jtau*P;
    
    ZZLB_DPE= 1/16*Rd*2*q(sqrt(M*SNR))+inv(J)*gammainc(M*SNR/2,3/2);
    ZZLB_tau=eye(M)* (1/8*varT*2*q(sqrt(SNR))+(1/(SNR*B_2*2))*gammainc(SNR/2,3/2));
    ZZLB_2SP=(P'*P)^-1*P'*(ZZLB_tau)*((P'*P)^-1*P')';
    
    
    fZZLB_DPE(SNRloop)=sqrt(ZZLB_DPE(1,1)+ZZLB_DPE(2,2)+ZZLB_DPE(3,3));
    fZZLB_2SP(SNRloop)=sqrt(ZZLB_2SP(1,1)+ZZLB_2SP(2,2)+ZZLB_2SP(3,3));
end