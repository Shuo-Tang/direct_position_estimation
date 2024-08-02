function [rawSignal_dpe] = readDataDPE_20ms(fid, numSol, settings)

if (settings.fileType==2) 
    dataAdaptCoeff=2;
end
rawSignal_dpe = zeros(numSol, 25375*20);
iMilliSec20 = 1;
iSol = 1;
hwb1 = waitbar(0,'Reading Raw Data...');
while iSol <= numSol
    % Read in the appropriate number of samples to process this
    % interation 
    [rawSignal, ~] = fread(fid, ...
        dataAdaptCoeff*25375*20, settings.dataType);

    rawSignal = rawSignal';
    rawSignal1=rawSignal(1:2:end);
    rawSignal2=rawSignal(2:2:end);
    rawSignal = rawSignal1 + 1i .* rawSignal2;  %transpose vector

    if mod(iMilliSec20 - 1, 50) == 0
        waitbar(iSol/numSol, hwb1, 'Reading Raw Data...');
        rawSignal_dpe(iSol,:) = rawSignal;
        iSol = iSol + 1;
    end
    iMilliSec20 = iMilliSec20 + 1; 
end
close(hwb1);
fclose(fid);
end

