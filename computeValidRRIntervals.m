function [RR_intervals, RR_times] = computeValidRRIntervals(heartlocs, invalidMask, fs)
% COMPUTEVALIDRRINTERVALS  RR intervals from R-peak sample indices, dropping
% any interval that spans an invalid (blanked) sample.
%
%   [RR_intervals, RR_times] = computeValidRRIntervals(heartlocs, invalidMask, fs)
%
% heartlocs   - R-peak sample indices (as saved in *_HRBR.mat)
% invalidMask - logical, sample-indexed, true = invalid/blanked (as saved
%               in *_HRBR.mat)
% fs          - sample rate (Hz)
%
% Extracted verbatim from HR_BR_HRVAnalysis_new.m's local function of the
% same name, so it can be reused directly against a *_HRBR.mat sidecar
% (heartlocs + invalidMask + fs, derivable from its saved t) without
% re-running the full raw-signal peak-detection pipeline -- the peak
% detection that produced heartlocs has already happened by the time an
% *_HRBR.mat file exists.

    if numel(heartlocs) < 2
        RR_intervals = []; RR_times = []; return;
    end
    csumInvalid  = [0; cumsum(invalidMask(:))];
    RR_intervals = [];
    RR_times     = [];
    for i = 1:numel(heartlocs)-1
        s1 = heartlocs(i);
        s2 = heartlocs(i+1);
        if s2 <= s1 + 1
            RR_intervals(end+1,1) = (s2 - s1) / fs; %#ok<AGROW>
            RR_times(end+1,1)     = s1 / fs;         %#ok<AGROW>
            continue;
        end
        if csumInvalid(s2) - csumInvalid(s1+1) == 0
            RR_intervals(end+1,1) = (s2 - s1) / fs; %#ok<AGROW>
            RR_times(end+1,1)     = s1 / fs;         %#ok<AGROW>
        end
    end
end
