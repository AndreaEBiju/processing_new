function [stimSegments, recoverySegments, stimOut, recoveryOut, stimIdx, recoveryIdx, stimHRChanIdx, recoveryHRChanIdx] = ...
    browseMotionArtifacts(data, fs, winSec, mode, baseFileName)
% browseMotionArtifacts
%
% Browse signals and enter motion artifact segments on close.
% Supports:
%   1) one-by-one dialogs
%   2) editable table
%   3) load from CSV
%   4) no artifacts option
%
% Also supports ABORTING execution cleanly.
%
% MODES
%   'single'         : data is [N x nChan]
%   'stim_recovery'  : data is a struct with fields .stim and .recovery

    if nargin < 3 || isempty(winSec),      winSec = 2;                    end
    if nargin < 4 || isempty(mode),        mode = 'single';               end
    if nargin < 5 || isempty(baseFileName),baseFileName = 'motion_artifacts'; end

    validateattributes(fs,     {'numeric'}, {'scalar','positive','real','finite'});
    validateattributes(winSec, {'numeric'}, {'scalar','positive','real','finite'});
    mode         = lower(string(mode));
    baseFileName = char(baseFileName);

    stimSegments     = []; recoverySegments = [];
    stimOut          = []; recoveryOut      = [];
    stimIdx          = []; recoveryIdx      = [];
    stimHRChanIdx    = []; recoveryHRChanIdx = [];

    switch mode
        case "single"
            validateattributes(data, {'numeric'}, {'2d','nonempty','real','finite'});
            [stimSegments, stimOut, stimIdx, stimHRChanIdx] = browseOneSignal( ...
                data, fs, winSec, 'Signal', ...
                [baseFileName '_segments.mat'], ...
                [baseFileName '_segment_indices.mat'], ...
                [baseFileName '_blankmotion.mat']);

        case "stim_rec"
            if ~isstruct(data) || ~isfield(data,'stim') || ~isfield(data,'recovery')
                error('For mode ''stim_recovery'', data must be a struct with fields .stim and .recovery');
            end
            validateattributes(data.stim,     {'numeric'}, {'2d','nonempty','real','finite'});
            validateattributes(data.recovery, {'numeric'}, {'2d','nonempty','real','finite'});

            [stimSegments, stimOut, stimIdx, stimHRChanIdx] = browseOneSignal( ...
                data.stim, fs, winSec, 'Stim Signal', ...
                [baseFileName '_stim_segments.mat'], ...
                [baseFileName '_stim_segment_indices.mat'], ...
                [baseFileName '_stim_blankmotion.mat']);

            if isempty(stimSegments) && isempty(stimIdx) && isempty(stimHRChanIdx)
                return;
            end

            [recoverySegments, recoveryOut, recoveryIdx, recoveryHRChanIdx] = browseOneSignal( ...
                data.recovery, fs, winSec, 'Recovery Signal', ...
                [baseFileName '_recovery_segments.mat'], ...
                [baseFileName '_recovery_segment_indices.mat'], ...
                [baseFileName '_recovery_blankmotion.mat']);

        otherwise
            error('Unknown mode: %s. Use ''single'' or ''stim_recovery''.', mode);
    end
end

% =========================================================================
function [removedSegments, yOut, removedSegmentIdx, hrChanIdx] = browseOneSignal( ...
    y, fs, winSec, figTitle, segmentFileName, idxFileName, blankmotionFileName)

    [N, nChan] = size(y);
    t          = (0:N-1) / fs;
    tEnd       = t(end);
    winDur     = winSec;                          % window width in seconds
    winSamp    = max(2, min(N, round(winSec*fs)));

    removedSegments    = [];
    removedSegmentIdx  = [];
    yOut               = y;
    hrChanIdx          = [];

    % ------------------------------------------------------------------
    % Adaptive display: ~5000 pts per channel in the visible window.
    % Full y is never modified.
    % ------------------------------------------------------------------
    N_DISPLAY = 5000;

    % ------------------------------------------------------------------
    % Build figure
    % ------------------------------------------------------------------
    fig = figure( ...
        'Name',            figTitle, ...
        'NumberTitle',     'off', ...
        'Units',           'normalized', ...
        'Position',        [0.08 0.08 0.84 0.84], ...
        'Color',           'w', ...
        'KeyPressFcn',     @keyScroll, ...
        'CloseRequestFcn', @onCloseRequest);

    ax    = gobjects(nChan, 1);
    hLine = gobjects(nChan, 1);

    [t_win, y_win] = getWindowData(0, winDur);

    for ch = 1:nChan
        ax(ch) = subplot(nChan, 1, ch, 'Parent', fig);
        hLine(ch) = plot(ax(ch), t_win, y_win(:,ch), 'k');
        ylabel(ax(ch), sprintf('Ch %d', ch));
        if ch < nChan
            ax(ch).XTickLabel = [];
        else
            xlabel(ax(ch), 'Time (s)');
        end
    end

    linkaxes(ax, 'x');
    xlim(ax(1), [0, min(tEnd, winDur)]);

    % ------------------------------------------------------------------
    % Slider — operates in seconds
    % ------------------------------------------------------------------
    sliderMax = max(winDur, tEnd - winDur);

    sld = uicontrol(fig, ...
        'Style',      'slider', ...
        'Units',      'normalized', ...
        'Position',   [0.10 0.02 0.80 0.03], ...
        'Min',        0, ...
        'Max',        sliderMax, ...
        'Value',      0, ...
        'SliderStep', [winDur/sliderMax, ...
                       winDur/sliderMax], ...
        'Callback',   @scrollCallback);

    uiwait(fig);

    % ------------------------------------------------------------------
    % Callbacks
    % ------------------------------------------------------------------
    function scrollCallback(src, ~)
        updatePlot(src.Value);
    end

    function keyScroll(~, event)
        step = winDur * 0.50;
        switch event.Key
            case 'rightarrow', sld.Value = min(sld.Value + step, sliderMax);
            case 'leftarrow',  sld.Value = max(sld.Value - step, 0);
            otherwise, return
        end
        updatePlot(sld.Value);
    end

    function updatePlot(t0)
        [t_win, y_win] = getWindowData(t0, winDur);
        for c = 1:nChan
            set(hLine(c), 'XData', t_win, 'YData', y_win(:,c));
        end
        xlim(ax(1), [t0, t0 + winDur]);
    end

    function [t_win, y_win] = getWindowData(t0, dur)
        i1  = max(1, floor(t0 * fs) + 1);
        i2  = min(N, ceil((t0 + dur) * fs) + 1);
        ds  = max(1, floor((i2 - i1 + 1) / N_DISPLAY));
        idx = i1:ds:i2;
        t_win = t(idx)';
        y_win = y(idx, :);
    end

    function onCloseRequest(~, ~)
        choice = questdlg( ...
            ['Before closing, note down the motion artifact time stamps.' newline newline ...
             'Choose how you want to proceed.'], ...
            'Close Motion Artifact Browser', ...
            'Continue Close', 'Keep Open', 'Abort', 'Keep Open');

        if isempty(choice) || strcmp(choice, 'Keep Open')
            return
        end

        if strcmp(choice, 'Abort')
            abortAndResume();
            return
        end

        % ---- ask about artifacts ----
        artifactChoice = questdlg( ...
            sprintf('For %s, do you want to mark artifact segments?', figTitle), ...
            'Artifact Choice', ...
            'Mark Segments', 'No Artifacts', 'Abort', 'Mark Segments');

        if isempty(artifactChoice), return, end

        if strcmp(artifactChoice, 'Abort')
            abortAndResume();
            return
        end

        if strcmp(artifactChoice, 'No Artifacts')
            removedSegments    = zeros(0,2);
            removedSegmentIdx  = zeros(0,2);
            yOut               = y;
            blankingApplied    = false;

            [hrChanIdx, hrStatus] = promptHRChannelIndex();
            if strcmp(hrStatus, 'abort')
                abortAndResume();
                return
            end
            if isempty(hrChanIdx), return, end

            save(segmentFileName,    'removedSegments');
            save(idxFileName,        'removedSegmentIdx');
            save(blankmotionFileName, ...
                'yOut','fs','t','removedSegments','removedSegmentIdx', ...
                'blankingApplied','hrChanIdx');
            closeAndResume();
            return
        end

        % ---- collect segments ----
        [segments, status] = collectSegmentsFlexible(tEnd, figTitle);

        if strcmp(status, 'abort')
            abortAndResume();
            return
        end
        if isempty(segments), return, end

        removedSegments    = segments;
        removedSegmentIdx  = timeToIndices(removedSegments, fs, N);

        save(segmentFileName, 'removedSegments');
        save(idxFileName,     'removedSegmentIdx');

        % ---- blanking ----
        blankChoice = questdlg( ...
            ['Apply motion-artifact blanking to ' figTitle '?' newline newline ...
             'Selected segments will be replaced with NaN.'], ...
            'Confirm Blanking', ...
            'Blank Signal', 'Do Not Blank', 'Abort', 'Do Not Blank');

        if isempty(blankChoice), return, end

        if strcmp(blankChoice, 'Abort')
            abortAndResume();
            return
        end

        blankingApplied = strcmp(blankChoice, 'Blank Signal');
        if blankingApplied
            yOut = applyBlanking(y, removedSegmentIdx);
        else
            yOut = y;
        end

        [hrChanIdx, hrStatus] = promptHRChannelIndex();
        if strcmp(hrStatus, 'abort')
            abortAndResume();
            return
        end
        if isempty(hrChanIdx), return, end

        save(blankmotionFileName, ...
            'yOut','fs','t','removedSegments','removedSegmentIdx', ...
            'blankingApplied','hrChanIdx');
        closeAndResume();
    end

    % ------------------------------------------------------------------
    % Helpers shared by onCloseRequest
    % ------------------------------------------------------------------
    function abortAndResume()
        removedSegments   = [];
        removedSegmentIdx = [];
        yOut              = y;
        hrChanIdx         = [];
        closeAndResume();
    end

    function closeAndResume()
        if isvalid(fig)
            uiresume(fig);
            delete(fig);
        end
    end
end

% =========================================================================
function [segments, status] = collectSegmentsFlexible(maxTime, signalLabel)
    segments = [];
    status   = 'keep_open';

    modeChoice = questdlg( ...
        sprintf('How do you want to enter artifact segments for %s?', signalLabel), ...
        'Artifact Segment Entry', ...
        'Edit Table', 'Load CSV', 'One by One', 'Edit Table');

    if isempty(modeChoice), return, end

    switch modeChoice
        case 'One by One', [segments, status] = collectSegmentsOneByOne(maxTime, signalLabel);
        case 'Edit Table',  [segments, status] = collectSegmentsTable(maxTime, signalLabel);
        case 'Load CSV',    [segments, status] = collectSegmentsCSV(maxTime, signalLabel);
    end
end

% =========================================================================
function [segments, status] = collectSegmentsOneByOne(maxTime, signalLabel)
    segments = [];
    status   = 'keep_open';

    while true
        answer = inputdlg( ...
            {'Start time (s):', 'Stop time (s):'}, ...
            sprintf('Enter segment to remove (%s)', signalLabel), ...
            [1 35; 1 35]);

        if isempty(answer)
            [segments, status] = abortOrKeepOpen('Cancel segment entry or abort the whole browse step?');
            return
        end

        tStart = str2double(strtrim(answer{1}));
        tStop  = str2double(strtrim(answer{2}));

        if ~validateSegmentMatrix([tStart, tStop], maxTime)
            uiwait(errordlg(sprintf('Times must satisfy:\n0 <= start < stop <= %.6g s', maxTime), ...
                'Invalid Segment', 'modal'));
            continue
        end

        segments = [segments; tStart, tStop]; %#ok<AGROW>

        nextChoice = questdlg( ...
            sprintf('Add another segment for %s?', signalLabel), ...
            'More Segments', ...
            'Add Another', 'Finish', 'Abort', 'Finish');

        if isempty(nextChoice)
            segments = []; status = 'keep_open'; return
        elseif strcmp(nextChoice, 'Abort')
            segments = []; status = 'abort';     return
        elseif strcmp(nextChoice, 'Finish')
            segments = sortrows(segments, 1);
            status   = 'ok';                     return
        end
    end
end

% =========================================================================
function [segments, status] = collectSegmentsCSV(maxTime, signalLabel)
    segments = [];
    status   = 'keep_open';

    [file, path] = uigetfile({'*.csv','CSV files (*.csv)'}, ...
        sprintf('Select CSV for %s', signalLabel));

    if isequal(file,0) || isequal(path,0)
        [~, status] = abortOrKeepOpen('CSV selection canceled. Keep browser open or abort?');
        return
    end

    T = readtable(fullfile(path, file));

    if width(T) < 2
        uiwait(errordlg('CSV must have at least 2 columns: start and stop.','Invalid CSV','modal'));
        return
    end

    segments = table2array(T(:,1:2));

    if ~validateSegmentMatrix(segments, maxTime)
        uiwait(errordlg(sprintf( ...
            'CSV contains invalid rows.\nEach row must satisfy: 0 <= start < stop <= %.6g s', maxTime), ...
            'Invalid CSV Segments', 'modal'));
        segments = []; return
    end

    segments = sortrows(segments, 1);
    status   = 'ok';
end

% =========================================================================
function [segments, status] = collectSegmentsTable(maxTime, signalLabel)
    segments = [];
    status   = 'keep_open';

    f = figure( ...
        'Name',        ['Artifact Segments: ' signalLabel], ...
        'NumberTitle', 'off', ...
        'MenuBar',     'none', ...
        'ToolBar',     'none', ...
        'Units',       'normalized', ...
        'Position',    [0.3 0.2 0.4 0.5], ...
        'WindowStyle', 'modal', ...
        'Color',       'w', ...
        'CloseRequestFcn', @closeTableFigure);

    uit = uitable(f, ...
        'Data',          nan(10,2), ...
        'ColumnName',    {'Start (s)', 'Stop (s)'}, ...
        'ColumnEditable',[true true], ...
        'Units',         'normalized', ...
        'Position',      [0.08 0.20 0.84 0.72]);

    uicontrol(f,'Style','pushbutton','String','Add 10 Rows', ...
        'Units','normalized','Position',[0.08 0.08 0.18 0.07], ...
        'Callback',@(~,~) set(uit,'Data',[uit.Data; nan(10,2)]));

    uicontrol(f,'Style','pushbutton','String','Load CSV', ...
        'Units','normalized','Position',[0.30 0.08 0.18 0.07], ...
        'Callback',@loadCSVCallback);

    uicontrol(f,'Style','pushbutton','String','OK', ...
        'Units','normalized','Position',[0.62 0.08 0.12 0.07], ...
        'Callback',@okCallback);

    uicontrol(f,'Style','pushbutton','String','Cancel', ...
        'Units','normalized','Position',[0.78 0.08 0.14 0.07], ...
        'Callback',@cancelCallback);

    uiwait(f);
    if isvalid(f), delete(f); end

    % --- nested callbacks ---
    function loadCSVCallback(~,~)
        [file, path] = uigetfile({'*.csv','CSV files (*.csv)'}, ...
            sprintf('Select CSV for %s', signalLabel));
        if isequal(file,0) || isequal(path,0), return, end
        T = readtable(fullfile(path, file));
        if width(T) < 2
            uiwait(errordlg('CSV must have at least 2 columns: start and stop.','Invalid CSV','modal'));
            return
        end
        uit.Data = table2array(T(:,1:2));
    end

    function okCallback(~,~)
        d    = uit.Data;
        keep = ~all(isnan(d), 2);
        d    = d(keep,:);

        if isempty(d)
            segments = []; status = 'keep_open';
            uiresume(f); return
        end

        if ~validateSegmentMatrix(d, maxTime)
            uiwait(errordlg(sprintf( ...
                'Invalid rows detected.\nEach row must satisfy: 0 <= start < stop <= %.6g s', maxTime), ...
                'Invalid Segments','modal'));
            return
        end

        segments = sortrows(d, 1);
        status   = 'ok';
        uiresume(f);
    end

    function cancelCallback(~,~)
        [segments, status] = abortOrKeepOpen('Cancel table entry or abort the whole browse step?');
        uiresume(f);
    end

    function closeTableFigure(~,~)
        cancelCallback();
    end
end

% =========================================================================
% Shared utilities
% =========================================================================
function [segments, status] = abortOrKeepOpen(msg)
    segments = [];
    choice = questdlg(msg, 'Cancel', 'Keep Open', 'Abort', 'Keep Open');
    if isempty(choice) || strcmp(choice, 'Keep Open')
        status = 'keep_open';
    else
        status = 'abort';
    end
end

function ok = validateSegmentMatrix(segments, maxTime)
    ok = true;
    if isempty(segments), return, end
    if ~isnumeric(segments) || size(segments,2) ~= 2, ok = false; return, end
    if any(isnan(segments(:))) || any(~isfinite(segments(:))), ok = false; return, end
    if any(segments(:,1) < 0) || any(segments(:,2) < 0) || ...
       any(segments(:,1) >= segments(:,2)) || any(segments(:,2) > maxTime)
        ok = false;
    end
end

function idxRanges = timeToIndices(timeRanges, fs, N)
    if isempty(timeRanges), idxRanges = zeros(0,2); return, end
    idxRanges = zeros(size(timeRanges));
    for k = 1:size(timeRanges,1)
        s = max(1, min(N, floor(timeRanges(k,1)*fs) + 1));
        e = max(1, min(N, ceil( timeRanges(k,2)*fs) + 1));
        if e < s, e = s; end
        idxRanges(k,:) = [s, e];
    end
end

function yBlanked = applyBlanking(y, idxRanges)
    yBlanked = y;
    for k = 1:size(idxRanges,1)
        yBlanked(idxRanges(k,1):idxRanges(k,2), :) = NaN;
    end
end

function [hrChanIdx, status] = promptHRChannelIndex()
    hrChanIdx = [];
    status    = 'keep_open';

    while true
        answer = inputdlg( ...
            {'Enter channel index for heart-rate detection (1-5):'}, ...
            'Heart Rate Channel', [1 40], {'1'});

        if isempty(answer)
            [~, status] = abortOrKeepOpen('Cancel HR channel entry or abort the whole browse step?');
            return
        end

        v = str2double(strtrim(answer{1}));

        if isnan(v) || ~isfinite(v) || v ~= round(v) || v < 1 || v > 5
            uiwait(errordlg('Channel index must be an integer between 1 and 5.', ...
                'Invalid Channel Index','modal'));
            continue
        end

        hrChanIdx = v;
        status    = 'ok';
        return
    end
end