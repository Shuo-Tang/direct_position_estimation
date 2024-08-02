function r = generateCAFs_baseline(config, rawSignal, ref_info, sat_info, caCode, pseudorange_est_ref, b)
% This function is designed for computing CAFs
% input:  usrInfo   -- receiver position, velocity, clock bias and clock
%                   drift information
%         satInfo   -- satellite position, velocity, clock bias, clock
%                   drift, el, az and PRN information
%         rawsignal -- received raw signal
%         caCode    -- C/A code generated for all GPS satellites based on
%                   sampling frequency
%         timeDiff  -- sampling moment difference w.r.t the minimum one 
%                      (if clock bias is searched, this will not be used)

%% reconstruct signal
x_ref_constructed = signalGenerate_given_baseline(config,...
                ref_info, sat_info, caCode, pseudorange_est_ref, b);

%% compute correlation
r = 0;
nSat = size(sat_info,1);
for iSat = 1:nSat
    r = r + abs(sum(x_ref_constructed(iSat,:).*conj(rawSignal)))^2;           
end
end




