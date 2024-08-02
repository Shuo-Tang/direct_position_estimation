function [rawSignal_dpe] = readDataDPE_1ms(fid, nSol, settings)
%% settins
if (settings.fileType == 8) 
    dataAdaptCoeff = 8;
end
samplesPerMs = settings.samplingFreq * 1e-3;
rawSignal_dpe = zeros(nSol, samplesPerMs);
iMilliSec = 1;
iSol = 1;
hwb1 = waitbar(0,'Reading Raw Data...');
while iSol <= nSol
    % Read the appropriate number of samples (1ms duration)
    [data, samplesRead] = fread(fid, [2, samplesPerMs], settings.dataType);
    rawSignal = data(1,:) + data(2,:).* 1i;

    if dataAdaptCoeff == 8 && (samplesRead ~= 2*samplesPerMs)
        disp('Not able to read the specified number of samples  for tracking, exiting!')
        fclose(fid);
        return
    end

    if mod(iMilliSec - 1, settings.navSolPeriod) == 0
        waitbar(iSol/nSol, hwb1, 'Reading Raw Data...');
        rawSignal_dpe(iSol,:) = rawSignal;
        iSol = iSol + 1;
    end
    iMilliSec = iMilliSec + 1; 
end
close(hwb1);
fclose(fid);
end

