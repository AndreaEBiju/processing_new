function outPath = pipeline_save_summary(D, P, conditionName)
% PIPELINE_SAVE_SUMMARY  Save a compact per-condition summary for comparison.
%
%   outPath = pipeline_save_summary(D, P, conditionName)
%
% Writes <neuralbase>_<condition>_summary.mat next to the source file. The
% summary holds, per neural channel, the binned noise-corrected excess RMS
% (step 3b) and the binned accepted-spike rate (steps 3-4) plus scalar means,
% the median noise floor, and the mean waveform. compare_conditions.m loads
% these across conditions.
%
% Prerequisites: steps 1 -> 4 (and 3b) have run, so D has .envelope and
% .spikes(k).alignedTimes.

    if nargin < 2 || isempty(P); P = pipeline_params(); end
    if nargin < 3 || isempty(conditionName); conditionName = 'cond'; end
    if ~isfield(D, 'envelope') || ~isfield(D, 'spikes')
        error('pipeline_save_summary:missing', 'Run steps through 3b/4 first.');
    end

    fs = D.fs; ch = D.neuralChannels; nCh = numel(ch);
    binSec = P.envBinSec;

    S = struct();
    S.condition = conditionName;
    S.fs        = fs;
    S.versionX  = D.versionX;
    S.labels    = D.channelLabels;
    S.channels  = ch;
    S.binSec    = binSec;
    S.wf_t_ms   = [];
    [S.t, S.excess_uv, S.rate_binned, S.validFrac, S.meanWaveform] = deal(cell(1, nCh));
    [S.featWidth, S.featVpp] = deal(cell(1, nCh));   % per-spike feature distributions
    [S.meanExcess, S.meanRate, S.medSigma, S.nSpikes, S.validSec] = deal(nan(1, nCh));
    [S.dfa_alpha1, S.dfa_alpha2, S.dfa_alphaFull, S.dfa_nCross] = deal(nan(1, nCh));
    S.dfa_nWindows = cell(1, nCh);
    hasDfa = isfield(D, 'dfa');

    for k = 1:nCh
        env = D.envelope(k);
        t   = env.t(:);
        S.t{k}          = t;
        S.excess_uv{k}  = env.excess_uv(:);
        S.validFrac{k}  = env.validFrac(:);
        S.meanExcess(k) = env.meanExcess_uv;

        % binned accepted-spike rate on the same bins as the envelope
        edges = [t - binSec/2; t(end) + binSec/2];
        st = D.spikes(k).alignedTimes(:);
        counts = histcounts(st, edges).';
        validSecBin = env.validFrac(:) * binSec;
        S.rate_binned{k} = counts ./ max(validSecBin, eps);

        S.meanRate(k)     = D.spikes(k).rate_hz;
        S.medSigma(k)     = median(D.sigma(D.validMask(:,k), k), 'omitnan') * 1e6;
        S.featWidth{k}    = D.spikes(k).width_ms(:);
        S.featVpp{k}      = D.spikes(k).Vpp_uv(:);
        S.meanWaveform{k} = D.spikes(k).meanWaveform;
        S.nSpikes(k)      = D.spikes(k).nSpikes;
        S.validSec(k)     = D.spikes(k).validSec;

        if hasDfa
            S.dfa_alpha1(k)    = D.dfa(k).alpha1;
            S.dfa_alpha2(k)    = D.dfa(k).alpha2;
            S.dfa_alphaFull(k) = D.dfa(k).alphaFull;
            S.dfa_nCross(k)    = D.dfa(k).nCross;
            S.dfa_nWindows{k}  = D.dfa(k).nWindows;
        end
    end
    if isfield(D.spikes, 'wf_t_ms'); S.wf_t_ms = D.spikes(1).wf_t_ms; end

    [folder, base] = fileparts(D.neuralFile);
    if isempty(folder); folder = pwd; end
    outPath = fullfile(folder, sprintf('%s_%s_summary.mat', base, conditionName));
    summary = S; %#ok<NASGU>
    save(outPath, 'summary');
    fprintf('[save] condition "%s" -> %s\n', conditionName, outPath);
    for k = 1:nCh
        fprintf('       ch %d (%s): mean excess %.2f uV | mean rate %.2f spk/s | sigma %.2f uV\n', ...
            ch(k), S.labels{k}, S.meanExcess(k), S.meanRate(k), S.medSigma(k));
    end
end
