function PosErrLS = conv2stepsPVT(r,sigen,config)

%% load configuration file
eval(config)

nfd = sigen.nfd;
fdCan = sigen.fdCan;
NsamplesLocal = sigen.NsamplesLocal;
UserPosition = param.UserPosition;
SatPosition = param.SatPosition;
deltaT = param.deltaT;
deltat = param.deltat;
%% memory allocation
EstRange = zeros(1,numSV);
% EstRxPVT = UserPosition'+10*(2*rand(3,1)-1);
% EstRxDeltaT = param.deltaT + 1e-8*rand;
x = zeros(4,num2stepsIterations);
x(:,1) = [UserPosition'+10*(2*rand(3,1)-1);(deltaT + 1e-8*rand)*c];
%% Estimate time delays and Doppler frequencies
fdEst = zeros(numSV,1);
delayEst = zeros(numSV,1);
RangeEst = zeros(numSV,1);
for kSV = 1:numSV
    [fdIndex,delayIndex] = find(r{kSV} == max(r{kSV}(:)));
    fdEst(kSV) = fdCan(fdIndex); % Doppler Frequency from Observation!
    delayEst(kSV) = (delayIndex-1)*dt;
    RangeEst(kSV) = delayEst(kSV)*c; % Range from Observation!
end
% visualization
% fdcan = (-5000:100:5000);
% delay = dt*(1:1:NsamplesLocal);
% figure, mesh(delay,fdcan,r{5}) 
%% Least Square Solution
for kIterations = 2:num2stepsIterations
    EstRxPVT = x(1:3,kIterations-1); 
    for kSV = 1:numSV
        EstRange(kSV) = norm(SatPosition(kSV,:)' - EstRxPVT(1:3)); % Range from Iterative Estimation!      
        numH = SatPosition(kSV, :)' - EstRxPVT(1:3);
        denH = norm(numH);
        H(kSV, 1:3) = -numH'/denH;
        H(kSV,4) = 1;   
    end   
    corrP  = (RangeEst - EstRange')/c;
    corrP_noAmbg = wrap(rem(corrP, 1e-3), 0.5e-3);   
    corrFracPseudorange = corrP_noAmbg*c;
    y = corrFracPseudorange + deltat'*c;
    correction = ((H'*H)\H')*y;
    x(1:4,kIterations) = x(1:4,kIterations-1) + correction;
end

PosErrLS=norm(x(1:3,end)-UserPosition');