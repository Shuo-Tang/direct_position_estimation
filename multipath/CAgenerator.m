function [ca_code]=CAgenerator(svnum)
k=length(svnum);            % number of satellites
ca_code=zeros(k,1023);
G1=zeros(1,1023);
G2=zeros(1,1023);
fs=4e6;
ts = 1/fs;   % Sampling period in sec
N=fs*1e-3;
numOfinterval=1000;
tc = 1/1.023e6;  % C/A chip period in sec
samplesPerCode=4000;
t=(1:N)/fs;
for i=1:k

  % the g2s vector holds the appropriate shift of the G2 code to
  % generate the C/A code (ex. for SV#19 -use a G2 shift of
  % g2s(19)=471) 
  
  g2s = [5;6;7;8;17;18;139;140;141;251;252;254;255;256;257;258;469;470;471;...
	 472;473;474;509;512;513;514;515;516; 859;860;861;862]; 
  g2shift=g2s(svnum(i),1); 
  % ***** Generate G1 code ***** % load shift register 
  reg= -1*ones(1,10); 
  for n= 1:1023
    G1(n) = reg(10) ; 
    save1 = reg(3)*reg(10); 			
    reg(1,2:10) = reg(1:1:9) ;
    reg(1) = save1; 
  end
  % ***** Generate G2 code ***** % load shift register 
  reg= -1*ones(1,10);
  for n= 1:1023
    g2(n)       = reg(10);
    saveBit     = reg(2)*reg(3)*reg(6)*reg(8)*reg(9)*reg(10);
    reg(2:10)   = reg(1:9);
    reg(1)      = saveBit;
  end 
  % ***** Shift G2 code to get G2i ***** 
  G2i =[g2(1023-g2shift+1 : 1023), g2(1 : 1023-g2shift)];
  % ***** Form single sample C/Acode by multiplying G1 and G2 
  ca_code= -(G1 .* G2i);
%   codeValueIndex = ceil((ts * (1:samplesPerCode)) / tc);
%   codeValueIndex(end) = 1023;
%   ca_code=ca_code(codeValueIndex);
end
end