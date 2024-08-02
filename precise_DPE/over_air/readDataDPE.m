function [rawSignal_dpe] = readDataDPE(fid, numSol, settings)

if (settings.fileType==2) 
    dataAdaptCoeff=2;
end
rawSignal_dpe = zeros(numSol, 25375);
iMilliSec = 1;
iSol = 1;
hwb1 = waitbar(0,'Reading Raw Data...');
while iSol <= numSol
    % Read in the appropriate number of samples to process this
    % interation 
    [rawSignal, ~] = fread(fid, ...
        dataAdaptCoeff*25375, settings.dataType);

    rawSignal = rawSignal';
    rawSignal1=rawSignal(1:2:end);
    rawSignal2=rawSignal(2:2:end);
    rawSignal = rawSignal1 + 1i .* rawSignal2;  %transpose vector

    if mod(iMilliSec - 1, 1000) == 0
        waitbar(iSol/numSol, hwb1, 'Reading Raw Data...');
        rawSignal_dpe(iSol,:) = rawSignal;
        iSol = iSol + 1;
    end
    iMilliSec = iMilliSec + 1; 
end
close(hwb1);
fclose(fid);
end

