naviBits = zeros(5, 2500*20);

for channelNr = activeChnList

    %=== Convert tracking output to navigation bits =======================

    %--- Copy 50000 long record from tracking output ---------------
    navBitsSamples = trackResults(channelNr).I_P';

    %--- Group every 20 vales of bits into columns ------------------------
    % navBitsSamples = reshape(navBitsSamples, ...
    %                          20, (size(navBitsSamples, 1) / 20));

    %--- Sum all samples in the bits to get the best estimate -------------
    % navBits = sum(navBitsSamples);

    %--- Now threshold and make 1 and 0 -----------------------------------
    % The expression (navBits > 0) returns an array with elements set to 1
    % if the condition is met and set to 0 if it is not met.
    navBits = (navBitsSamples > 0);

    %--- Convert from decimal to binary -----------------------------------
    % The function ephemeris expects input in binary form. In Matlab it is
    % a string array containing only "0" and "1" characters.
    navBitsBin = dec2bin(navBits);
    
    %--- Record the navi bits -----------------------------------
    naviBits(channelNr, :) = navBits;
end

% naviBitsExtended = kron(naviBits, ones(1,20));
% naviBitsExtended(naviBitsExtended == 0) = -1;
% save("naviBitsExtended.mat", "naviBitsExtended")

naviBits(naviBits == 0) = -1;
save("naviBits.mat", "naviBits")