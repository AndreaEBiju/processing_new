function P = pipeline_params()
% PIPELINE_PARAMS  Central tunable parameters for the VENG analysis pipeline.
%
% Every step reads its knobs from this struct so that the same settings are
% used in single-file (interactive) and bulk (batch) modes. Edit values here
% rather than inside the step functions.

    P = struct();

    % ---- Dataset / naming (used by the bulk loader later) -----------------
    % Files are expected to carry a version tag like "v0.x.y" in their name.
    % In bulk mode only files whose middle number matches P.versionX are
    % processed. Leave empty to accept any version in single-file mode.
    P.versionX           = [];     % e.g. 3  -> matches v0.3.*  (set in GUI/bulk)
    P.useSlowWave        = false;  % ask for / load the slow-wave file in the loader
                                   % (only needed for slow-wave phase-locking; off for now)

    % ---- Step 1: bandpass filter (the spike band) -------------------------
    P.bandpassLow        = 100;    % Hz  high-pass corner (removes ECG/resp drift)
    P.bandpassHigh       = 5000;   % Hz  low-pass corner (clamped to fs/2-1)
    P.filterOrder        = 4;      % Butterworth order (zero-phase via filtfilt)

    % ---- Step 1b: cardiac (QRS) template subtraction ---------------------
    P.cardiacPreMs        = 10;    % template window before R-fiducial (ms)
    P.cardiacPostMs       = 20;    % template window after R-fiducial (ms)
    P.cardiacSearchMs     = 5;     % search +/- this around each R-peak for the QRS peak
    P.cardiacMinActiveFrac = 0.05; % skip beats whose window amplitude < this x median
                                   % (keeps the dead-zone from being touched)
    P.cardiacTemplateBeats = 40;   % rolling template = median of nearest K beats
                                   % (tracks slow QRS-shape drift; large K -> global)

    % ---- Step 2: noise estimate + validity mask --------------------------
    P.sigmaWindowSec     = 5;      % sliding window length for robust noise sigma
    P.sigmaStepFrac      = 0.5;    % window step as a fraction of the window (overlap)
    P.edgeBufferMs       = 5;      % pad (ms) around invalid regions (filter ringing)
    P.zeroRunMinSec      = 0.5;    % min duration (s) of a flat/dead run to flag invalid
    P.flatRangeFrac      = 0.05;   % dead if local range < this fraction of the
                                   % recording's typical local range (catches
                                   % near-zero dropouts, not just exact zeros)
    P.sigmaMinValidFrac  = 0.2;    % min valid fraction in a window to estimate sigma

    % ---- Step 3b: spike-band activity envelope (RMS) ---------------------
    P.envBinSec          = 1;      % RMS bin width (s) for the activity trace
    P.envCardiacGuardMs  = 15;     % exclude +/- this around each R-peak from RMS
                                   % (so residual cardiac doesn't inflate activity)
    P.envSmoothSec       = 0.5;    % smoothing for the continuous display envelope
    P.envArtifactSigma   = 15;     % exclude excursions above this x sigma from the RMS
                                   % (unguarded cardiac / motion; neural is < ~3 sigma)

    % ---- Step 3: spike detection -----------------------------------------
    P.threshSigma        = 4.5;    % detection threshold in units of sigma
    P.detectPolarity     = 'neg';  % 'neg' | 'pos' | 'both'
    P.refractoryMs       = 1.0;    % dead time between detections (ms)
    P.maxThreshSigma     = 40;     % peaks above this are flagged as artifacts

    % ---- Step 4: waveform extraction + alignment -------------------------
    P.wfPreMs            = 1.0;    % waveform window before the aligned peak (ms)
    P.wfPostMs           = 2.0;    % waveform window after the aligned peak (ms)
    P.wfAlignSearchMs    = 0.5;    % re-align to the extremum within +/- this (ms)
    % waveform screening (reject noise crossings / artifact residual)
    P.minAmpUV           = 8;      % min peak-to-peak (uV)
    P.maxAmpUV           = 150;    % max Vpp; above this = cardiac/motion artifact
    P.minWidthMs         = 0.2;    % min FWHM; below = single-sample noise
    P.maxWidthMs         = 2.5;    % max FWHM; above = broadband / window-edge noise

    % ---- Step 5: sort check (feature space + density clustering) ---------
    P.numPCs             = 6;      % PCA components kept as clustering features
    P.tsnePerplexity     = 30;     % t-SNE perplexity (auto-capped to n/4)
    P.dbscanMinPts       = 10;     % DBSCAN min points per cluster
    P.dbscanEps          = [];     % DBSCAN epsilon ([] -> auto from k-distance)
    P.dbscanEpsPct       = 75;     % percentile of k-distance used when eps is auto

    % ---- Step 9: cardiac-locked spike removal (peri-R-wave) --------------
    P.cardiacRemove      = true;   % remove cardiac-locked spikes + censor their time
    P.cardiacUniform     = true;   % apply a FIXED window to ALL channels (uniform), not
                                   % only coupled ones -> comparable across channels/conditions
    P.cardiacRemoveWinMs = 15;     % PRIMARY USE: half-window (ms) NaN-blanked around each
                                   % R-peak in step1a (blank-before-filter). Also reused by the
                                   % step9 tag-only QC diagnostic as its exclusion half-window.
    P.cardiacLagMaxMs    = 75;     % peri-R histogram half-range (ms)
    P.cardiacLagBinMs    = 2;      % peri-R histogram bin (ms)
    P.cardiacBaseLagMs   = 40;     % baseline from |lag| > this (ms)
    P.cardiacZThresh     = 4;      % min peak z (vs baseline) to call coupling
    P.cardiacWinMaxMs    = 30;     % cap on the exclusion-window half-extent (ms)

    % ---- Step 5b: population verification (over-cluster + d' merge) -------
    P.verifyKover         = 12;    % initial over-clustering k (k-means)
    P.verifyMinDprime     = 2.0;   % merge populations whose templates are separated
                                   % by less than this d' (i.e. overlap within noise)
    P.verifyMinClusterSize = 20;   % drop populations smaller than this after merging

    % ---- Step 5c: modality test (Silverman smoothed-bootstrap) -----------
    P.modalityNBoot      = 200;    % bootstrap replicates for the unimodality p-value
    P.modalityGrid       = 512;    % KDE grid points
    P.modalityMaxN       = 2000;   % subsample to this many spikes for speed
    P.modalityAlpha      = 0.05;   % significance for declaring multimodal

    % ---- Step 6: single-population summary / spike-train report ----------
    P.frBinSec           = 1;      % firing-rate bin width (s)
    P.rasterRowSec       = 60;     % raster wrap: one row per this many seconds
    P.wfBandPct          = [10 90];% waveform shaded band percentiles (range)
    P.vppBinSec          = 5;      % bin for peak-to-peak amplitude over time
    P.acgMaxLagSec       = 0.5;    % autocorrelogram lag range (+/-)
    P.acgBinSec          = 0.005;  % autocorrelogram bin
    P.fanoMinWinSec      = 0.02;   % Fano: smallest counting window
    P.fanoMaxWinSec      = 5;      % Fano: largest counting window
    P.fanoNWin           = 15;     % Fano: number of window sizes (log-spaced)
    P.fanoCanonSec       = 0.1;    % canonical Fano window (s) reported as a scalar
                                   % (short, so fully-valid windows exist under censoring)
    P.psdRateBinSec      = 0.05;   % firing-rate sampling for power spectrum (20 Hz)
    P.psdMaxHz           = 10;     % display up to this frequency
    P.cv2WinSec          = 30;     % window for rolling CV2 over time
    P.isiHeatTBinSec     = 10;     % time bin for the time x ISI heatmap
    % burst detection (data-driven log-ISI / void parameter)
    P.burstMinSpikes     = 3;      % min spikes per burst
    P.burstVoidThresh    = 0.7;    % bimodality required (Pasquale void parameter)
    P.burstMaxThreshMs   = 100;    % intra-burst ISI peak must be below this

    % ---- Plotting ---------------------------------------------------------
    P.plotMaxPoints      = 200000; % max points drawn in full-trace plots
                                   % (display decimation only; analysis uses
                                   %  the full-resolution signal)

end
