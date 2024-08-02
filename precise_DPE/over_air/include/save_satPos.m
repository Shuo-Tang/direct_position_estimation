
nEpoch = length(SatelliteInformation);
PRNlist = [10, 27, 8, 16];
nSat = length(PRNlist);

satPos_navFile = zeros(nSat, 3, nEpoch);

for iEpoch = 1:nEpoch
    allSat = SatelliteInformation(iEpoch).labelSatRover{1, 1};
    allSatPos = SatelliteInformation(iEpoch).posSatRover{1, 1};
    for iSat = 1:nSat
        satPRN = PRNlist(iSat);
        satIndex = find(allSat == satPRN);
        satPos_navFile(iSat, :, iEpoch) = allSatPos(satIndex,:);
    end
end

save("satPos_navFile.mat", "satPos_navFile")