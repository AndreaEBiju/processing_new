function [y, blankIdx, artifact_locs] = blank_stim_spikes_nan(signal, fs)
% blank_stim_spikes_nan_allch
% Detect stimulation artifacts across all channels, merge detections,
% and blank the same time windows across all channels using NaN.
%
% INPUTS
%   signal : [N x nChan] signal matrix
%   fs     : sampling rate
%
% OUTPUTS
%   y             : signal with stim artifact windows set to NaN
%   blankIdx      : [nArtifacts x 2] sample ranges blanked
%   artifact_locs : detected artifact center samples

    pre_blank_ms = 5;
    post_blank_ms = 5;
    artifact_min_distance_ms = 5;

    pre_blank_samples = round(pre_blank_ms/1000 * fs);
    post_blank_samples = round(post_blank_ms/1000 * fs);
    artifact_min_distance_samples = round(artifact_min_distance_ms/1000 * fs);

    [N, nChan] = size(signal);

    all_locs = [];

    for ch = 1:nChan
        x = signal(:, ch);

        if all(isnan(x)) || all(x == 0 | isnan(x))
            continue;
        end

        artifact_threshold = 0.5 * max(abs(x), [], 'omitnan');

        if ~isfinite(artifact_threshold) || artifact_threshold <= 0
            continue;
        end

        [~, locs] = findpeaks(abs(x), ...
            'MinPeakHeight', artifact_threshold, ...
            'MinPeakDistance', artifact_min_distance_samples);

        all_locs = [all_locs; locs(:)]; %#ok<AGROW>
    end

    if isempty(all_locs)
        y = signal;
        blankIdx = zeros(0,2);
        artifact_locs = [];
        return;
    end

    all_locs = sort(all_locs);

    % Merge artifact detections across channels
    artifact_locs = all_locs(1);

    for i = 2:numel(all_locs)
        if all_locs(i) - artifact_locs(end) > artifact_min_distance_samples
            artifact_locs(end+1,1) = all_locs(i); %#ok<AGROW>
        else
            artifact_locs(end) = round(mean([artifact_locs(end), all_locs(i)]));
        end
    end

    y = signal;
    blankIdx = zeros(numel(artifact_locs), 2);

    for k = 1:numel(artifact_locs)
        s1 = max(1, artifact_locs(k) - pre_blank_samples);
        s2 = min(N, artifact_locs(k) + post_blank_samples);

        blankIdx(k,:) = [s1 s2];
        y(s1:s2, :) = NaN;
    end
end