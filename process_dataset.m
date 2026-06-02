function D = process_dataset(D, P, plotMode)
% PROCESS_DATASET  Run the full per-file pipeline (steps 1 -> 5c) on one D.
%
%   D = process_dataset(D, P, plotMode)
%
% Reusable so the same sequence runs in single-file, pair, and bulk modes.
% plotMode = true draws every step's diagnostic figures; false (bulk/pair)
% suppresses them. Steps 5 and 5b (exploratory cluster views) are NOT run
% here -- call them manually if wanted; 5c (modality verdict) and 6 (the
% spike-train report) are included. step6 always computes D.metrics; it only
% draws the 12-panel report figures when plotMode is true.

    if nargin < 2 || isempty(P);        P = pipeline_params(); end
    if nargin < 3 || isempty(plotMode); plotMode = true;       end

    D = step1a_blank_cardiac(D, P, plotMode);   % NaN-blank R-peak windows pre-filter
    D = step1_bandpass(D, P, plotMode);         % interpolates across blanks -> no ringing
    D = step2_noise_sigma(D, P, plotMode);
    D = step3_detect(D, P, plotMode);
    D = step3b_envelope(D, P, plotMode);
    D = step4_waveforms(D, P, plotMode);
    D = step5c_modality_test(D, P, plotMode);
    D = step6_spike_report(D, P, plotMode);
    % NOTE: cardiac is handled by blanking before the filter (step1a); the old
    % template subtraction (step1b) and post-hoc removal (step9) are no longer
    % in the flow. step9 remains available as a tag-only diagnostic for the
    % peri-R contamination figure (run it on an un-blanked detection).
end
