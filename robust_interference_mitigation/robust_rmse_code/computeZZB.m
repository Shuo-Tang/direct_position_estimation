x_local = sigen.x_local;
% iSRN = 7;
% B_2=sum((diff(x_local(iSRN,:))/dt).^2)/sum(x_local(iSRN,:).^2);
B_2=sum((diff(x_local')'/dt).^2, 2)./sum(x_local.^2, 2);
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
    if RIMuse
        Th_norm = 1e-15;
        LoE = (1 - exp(-Th_norm^2 + Th_norm*sqrt(pi)/2*erfc(Th_norm)))^2/(1 - exp(-Th_norm^2));
        SNR = SNR*LoE;
    end
    Jtau=diag(repmat(SNR*2,1,M).*B_2');
    J= P'*Jtau*P;
    B_2 = B_2(1);
    ZZLB_DPE= 1/16*Rd*2*q(sqrt(M*SNR))+inv(J)*gammainc(M*SNR/2,3/2);
    ZZLB_tau=eye(M)* (1/8*varT*2*q(sqrt(SNR))+(1/(SNR*B_2*2))*gammainc(SNR/2,3/2));
    ZZLB_2SP=(P'*P)^-1*P'*(ZZLB_tau)*((P'*P)^-1*P')';
    
    
    fZZLB_DPE(SNRloop)=sqrt(ZZLB_DPE(1,1)+ZZLB_DPE(2,2)+ZZLB_DPE(3,3));
    fZZLB_2SP(SNRloop)=sqrt(ZZLB_2SP(1,1)+ZZLB_2SP(2,2)+ZZLB_2SP(3,3));
end