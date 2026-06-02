function D = step0_load_data(P, titleSuffix)
% STEP0_LOAD_DATA  Interactive loader for one VENG dataset (single-file mode).
%
%   D = step0_load_data(P)
%
% Opens a GUI to:
%   * browse and select the neural data .mat file;
%   * map the variables inside it (defaults to the earlier schema:
%     yOut / fs / t / removedSegmentIdx, all overridable);
%   * mark each column of the neural matrix as Neural / SlowWaveRaw / Ignore
%     and give the neural channels labels;
%   * browse and select the heartbeat file (R-peaks read as SAMPLE INDICES);
%   * browse and select the extracted slow-wave file (3 channels at neural fs);
%   * set the version number x (informational here; used for matching in bulk).
%
% Returns a struct D consumed by the later step functions, or [] if cancelled.
%
% Notes
% -----
% * Stim/recovery splitting is assumed ALREADY DONE upstream; this loader
%   does not split.
% * Raw slow-wave columns in the neural matrix are NOT used for slow-wave
%   analysis (mark them Ignore); slow-wave data comes from the separate
%   extracted file, which is expected to be 3 channels sampled at neural fs.

    if nargin < 1 || isempty(P); P = pipeline_params(); end
    if nargin < 2 || isempty(titleSuffix); titleSuffix = ''; end

    D = [];   % default = cancelled

    % ---- shared state held in this function's workspace -------------------
    st = struct();
    st.neuralPath = '';  st.neuralVars = {};  st.neuralSizes = {};
    st.rpeakPath  = '';  st.rpeakVars  = {};
    st.swPath     = '';  st.swVars     = {};
    st.cancelled  = false;

    % ====================== build figure ===================================
    fig = uifigure('Name', ['Step 0 — Load dataset' titleSuffix], 'Position', [180 120 760 720]);
    g = uigridlayout(fig, [1 1]); g.Padding = [10 10 10 10];
    sp = uipanel(g, 'BorderType', 'none', 'Scrollable', 'on');
    L = uigridlayout(sp, [20 4]);
    L.RowHeight   = repmat({'fit'}, 1, 20);
    L.ColumnWidth = {150, '1x', 150, '1x'};
    L.RowSpacing  = 6; L.ColumnSpacing = 8; L.Padding = [6 6 6 6];

    r = 1;
    sectionLabel(L, r, 'NEURAL DATA FILE'); r = r + 1;

    btnNeural = uibutton(L, 'Text', 'Browse neural .mat...');
    btnNeural.Layout.Row = r; btnNeural.Layout.Column = 1;
    lblNeural = uilabel(L, 'Text', '(none selected)', 'FontColor', [0.4 0.4 0.4]);
    lblNeural.Layout.Row = r; lblNeural.Layout.Column = [2 4]; r = r + 1;

    [ddNeuralMat, r] = mapRow(L, r, 'Neural matrix var:', 'fs var:');
    ddNeuralMatVar = ddNeuralMat{1}; ddFsVar = ddNeuralMat{2};
    [ddNeural2, r]  = mapRow(L, r, 'Time (t) var:', 'removedSegmentIdx var:');
    ddTVar = ddNeural2{1}; ddSegVar = ddNeural2{2};

    chanHdr = uilabel(L, 'Text', 'Channel assignment (mark raw slow-wave columns as Ignore):', ...
        'FontAngle', 'italic');
    chanHdr.Layout.Row = r; chanHdr.Layout.Column = [1 4]; r = r + 1;

    chanTbl = uitable(L);
    chanTbl.Layout.Row = r; chanTbl.Layout.Column = [1 4];
    L.RowHeight{r} = 150;
    chanTbl.ColumnName     = {'Column', 'Type', 'Label'};
    chanTbl.ColumnEditable = [false true true];
    chanTbl.Data = emptyChanTable();
    r = r + 1;

    sectionLabel(L, r, 'HEARTBEAT (R-PEAK) FILE  —  read as sample indices'); r = r + 1;
    btnR = uibutton(L, 'Text', 'Browse R-peak .mat...');
    btnR.Layout.Row = r; btnR.Layout.Column = 1;
    lblR = uilabel(L, 'Text', '(none selected)', 'FontColor', [0.4 0.4 0.4]);
    lblR.Layout.Row = r; lblR.Layout.Column = [2 4]; r = r + 1;
    [ddRcell, r] = mapRow(L, r, 'R-peak var:', '');
    ddRVar = ddRcell{1};

    btnSW = []; lblSW = []; ddSWVar = [];
    if P.useSlowWave
        sectionLabel(L, r, 'EXTRACTED SLOW-WAVE FILE  —  3 channels @ neural fs'); r = r + 1;
        btnSW = uibutton(L, 'Text', 'Browse slow-wave .mat...');
        btnSW.Layout.Row = r; btnSW.Layout.Column = 1;
        lblSW = uilabel(L, 'Text', '(none selected)', 'FontColor', [0.4 0.4 0.4]);
        lblSW.Layout.Row = r; lblSW.Layout.Column = [2 4]; r = r + 1;
        [ddSWcell, r] = mapRow(L, r, 'Slow-wave var:', '');
        ddSWVar = ddSWcell{1};
    end

    sectionLabel(L, r, 'DATASET'); r = r + 1;
    vlbl = uilabel(L, 'Text', 'Version x (v0.x.y):');
    vlbl.Layout.Row = r; vlbl.Layout.Column = 1;
    efVersion = uieditfield(L, 'numeric', 'Value', 0, 'Limits', [0 Inf], ...
        'RoundFractionalValues', true);
    efVersion.Layout.Row = r; efVersion.Layout.Column = 2;
    if ~isempty(P.versionX); efVersion.Value = P.versionX; end
    r = r + 1;

    btnGo = uibutton(L, 'Text', 'Continue', 'BackgroundColor', [0.7 0.9 0.7]);
    btnGo.Layout.Row = r; btnGo.Layout.Column = 4;
    status = uilabel(L, 'Text', '', 'FontColor', [0.7 0 0]);
    status.Layout.Row = r; status.Layout.Column = [1 3];

    % ---- callbacks --------------------------------------------------------
    btnNeural.ButtonPushedFcn = @(s,e) onBrowseNeural();
    btnR.ButtonPushedFcn      = @(s,e) onBrowse('rpeak', lblR, {ddRVar});
    if P.useSlowWave
        btnSW.ButtonPushedFcn = @(s,e) onBrowse('sw', lblSW, {ddSWVar});
    end
    ddNeuralMatVar.ValueChangedFcn = @(s,e) refreshChannelTable();
    btnGo.ButtonPushedFcn     = @(s,e) onContinue();
    fig.CloseRequestFcn       = @(s,e) onCancel();

    uiwait(fig);
    if isvalid(fig); delete(fig); end
    return;

    % ====================== nested helpers =================================
    function onBrowseNeural()
        [f, p] = uigetfile({'*.mat','MAT-files'}, 'Select neural data file', pref_dir('neural'));
        raiseFig();                        % bring the loader back to the front
        if isequal(f,0); return; end
        full = fullfile(p, f);
        lblNeural.Text = ['reading ' f ' ...']; status.Text = ''; drawnow;
        try
            info = whos('-file', full);
        catch ME
            lblNeural.Text = '(none selected)';
            status.Text = ['Cannot read "' f '": ' ME.message];
            return;
        end
        if isempty(info)
            lblNeural.Text = '(none selected)';
            status.Text = ['No variables found in "' f '".']; return;
        end
        save_pref_dir('neural', p);
        st.neuralPath = full;
        lblNeural.Text = f;
        st.neuralVars  = {info.name};
        st.neuralSizes = {info.size};
        setDropdown(ddNeuralMatVar, st.neuralVars, 'yOut');
        setDropdown(ddFsVar,        st.neuralVars, 'fs');
        setDropdown(ddTVar,    [{'<none>'} st.neuralVars], 't');
        setDropdown(ddSegVar,  [{'<none>'} st.neuralVars], 'removedSegmentIdx');
        refreshChannelTable();
    end

    function onBrowse(kind, lbl, dds)
        [f, p] = uigetfile({'*.mat','MAT-files'}, 'Select file', pref_dir(kind));
        raiseFig();                        % bring the loader back to the front
        if isequal(f,0); return; end
        full = fullfile(p, f);
        % feedback: reading the header can block (esp. cloud/large files)
        lbl.Text = ['reading ' f ' ...']; status.Text = ''; drawnow;
        try
            info = whos('-file', full);
        catch ME
            lbl.Text = '(none selected)';
            status.Text = ['Cannot read "' f '": ' ME.message];
            return;
        end
        if isempty(info)
            lbl.Text = '(none selected)';
            status.Text = ['No variables found in "' f '".'];
            return;
        end
        save_pref_dir(kind, p);
        vars = {info.name};
        lbl.Text = f;
        switch kind
            case 'rpeak'
                st.rpeakPath = full; st.rpeakVars = vars;
                setDropdown(dds{1}, vars, guessVector(info));
            case 'sw'
                st.swPath = full; st.swVars = vars;
                setDropdown(dds{1}, vars, guessMatrix(info));
        end
    end

    function refreshChannelTable()
        if isempty(st.neuralPath); return; end
        v = ddNeuralMatVar.Value;
        idx = find(strcmp(st.neuralVars, v), 1);
        if isempty(idx); chanTbl.Data = emptyChanTable(); return; end
        sz = st.neuralSizes{idx};
        % columns = channels (samples x channels). Use the smaller dim as #ch.
        if numel(sz) < 2; chanTbl.Data = emptyChanTable(); return; end
        nCh = min(sz);   % robust to orientation
        col = (1:nCh)';
        type = categorical(repmat({'Neural'}, nCh, 1), {'Neural','SlowWaveRaw','Ignore'});
        lab  = strings(nCh, 1);
        for i = 1:nCh; lab(i) = sprintf('Ch%d', i); end
        chanTbl.Data = table(col, type, lab, 'VariableNames', {'Column','Type','Label'});
    end

    function onContinue()
        status.Text = '';
        if isempty(st.neuralPath); status.Text = 'Select a neural file.'; return; end
        if isempty(st.rpeakPath);  status.Text = 'Select an R-peak file.'; return; end
        if P.useSlowWave && isempty(st.swPath); status.Text = 'Select a slow-wave file.'; return; end
        btnGo.Enable = 'off'; btnGo.Text = 'Loading...';
        dlg = uiprogressdlg(fig, 'Title', 'Step 0', ...
            'Message', 'Loading files (a large neural file may take a while)...', ...
            'Indeterminate', 'on');
        drawnow;
        try
            D = doLoad();
        catch ME
            close(dlg);
            btnGo.Enable = 'on'; btnGo.Text = 'Continue';
            status.Text = ['Load error: ' ME.message];
            return;
        end
        close(dlg);
        uiresume(fig);
    end

    function onCancel()
        st.cancelled = true;
        D = [];
        uiresume(fig);
    end

    function raiseFig()
        % keep the loader on top after a uigetfile dialog closes
        try
            figure(fig);
        catch
            try; fig.Visible = 'off'; fig.Visible = 'on'; catch; end
        end
    end

    % ====================== the actual load ================================
    function Dout = doLoad()
        yVar  = ddNeuralMatVar.Value;
        fsVar = ddFsVar.Value;

        % Load only the variables we actually need. Loading the whole file is
        % what made Continue feel frozen when the file holds extra big arrays.
        nVars = {yVar, fsVar};
        if ~strcmp(ddTVar.Value, '<none>');   nVars{end+1} = ddTVar.Value;   end
        if ~strcmp(ddSegVar.Value, '<none>'); nVars{end+1} = ddSegVar.Value; end
        fprintf('[step0] Loading neural variables: %s ...\n', strjoin(nVars, ', '));
        Sn = load(st.neuralPath, nVars{:});
        if ~isfield(Sn, yVar);  error('Variable "%s" not found in neural file.', yVar);  end
        if ~isfield(Sn, fsVar); error('Variable "%s" not found in neural file.', fsVar); end

        y  = double(Sn.(yVar));
        fs = double(Sn.(fsVar));
        fs = fs(1);

        % orient neural matrix to samples x channels (more rows than cols)
        if size(y,1) < size(y,2); y = y.'; end

        % optional time vector
        if ~strcmp(ddTVar.Value, '<none>') && isfield(Sn, ddTVar.Value)
            t = Sn.(ddTVar.Value); t = t(:);
        else
            t = (0:size(y,1)-1)' / fs;
        end

        % optional blanked-segment indices
        if ~strcmp(ddSegVar.Value, '<none>') && isfield(Sn, ddSegVar.Value)
            removedSegmentIdx = Sn.(ddSegVar.Value);
        else
            removedSegmentIdx = zeros(0, 2);
        end

        % channel assignment
        Tc = chanTbl.Data;
        types = string(Tc.Type);
        neuralChannels = Tc.Column(types == "Neural")';
        labels = cellstr(Tc.Label(types == "Neural"));
        if isempty(neuralChannels)
            error('No channels marked as Neural.');
        end
        if max(neuralChannels) > size(y,2)
            error('Neural channel index exceeds matrix columns.');
        end

        % R-peaks as SAMPLE INDICES
        Sr = load(st.rpeakPath, ddRVar.Value);
        if ~isfield(Sr, ddRVar.Value)
            error('R-peak variable "%s" not found.', ddRVar.Value);
        end
        rpeakSamples = round(double(Sr.(ddRVar.Value)(:)));
        rpeakSamples = rpeakSamples(rpeakSamples >= 1 & rpeakSamples <= size(y,1));
        rpeakTimes   = (rpeakSamples - 1) / fs;

        % extracted slow-wave (optional): 3 channels at neural fs
        sw = [];
        if P.useSlowWave
            Sw = load(st.swPath, ddSWVar.Value);
            if ~isfield(Sw, ddSWVar.Value)
                error('Slow-wave variable "%s" not found.', ddSWVar.Value);
            end
            sw = double(Sw.(ddSWVar.Value));
            if size(sw,1) < size(sw,2); sw = sw.'; end   % samples x channels
            if size(sw,1) ~= size(y,1)
                warning('step0:swLength', ...
                    'Slow-wave length (%d) ~= neural length (%d); using as-is.', ...
                    size(sw,1), size(y,1));
            end
        end

        Dout = struct();
        Dout.neuralFile        = st.neuralPath;
        Dout.rpeakFile         = st.rpeakPath;
        Dout.slowWaveFile      = st.swPath;
        Dout.versionX          = efVersion.Value;
        Dout.fs                = fs;
        Dout.y                 = y;                 % samples x channels (raw)
        Dout.t                 = t;
        Dout.removedSegmentIdx = removedSegmentIdx;
        Dout.neuralChannels    = neuralChannels;    % column indices into y
        Dout.channelLabels     = labels;
        Dout.rpeakSamples      = rpeakSamples;
        Dout.rpeakTimes        = rpeakTimes;
        Dout.slowWave          = sw;                % samples x 3 @ neural fs
        Dout.meta              = struct('loadedAt', datetime('now'), ...
                                        'neuralVar', yVar, 'fsVar', fsVar);

        fprintf('[step0] Loaded %s\n', st.neuralPath);
        fprintf('        fs = %g Hz | %d samples (%.1f s) | %d neural channel(s): %s\n', ...
            fs, size(y,1), size(y,1)/fs, numel(neuralChannels), strjoin(labels, ', '));
        if P.useSlowWave
            fprintf('        %d R-peaks | slow-wave %d x %d\n', ...
                numel(rpeakSamples), size(sw,1), size(sw,2));
        else
            fprintf('        %d R-peaks | slow-wave: not loaded\n', numel(rpeakSamples));
        end
    end
end

% ========================= file-scope helpers =============================
function sectionLabel(L, row, txt)
    lb = uilabel(L, 'Text', txt, 'FontWeight', 'bold', 'FontColor', [0.15 0.30 0.55]);
    lb.Layout.Row = row; lb.Layout.Column = [1 4];
end

function [dds, rowOut] = mapRow(L, row, lab1, lab2)
    u1 = uilabel(L, 'Text', lab1); u1.Layout.Row = row; u1.Layout.Column = 1;
    d1 = uidropdown(L, 'Items', {'<none>'}, 'Value', '<none>');
    d1.Layout.Row = row; d1.Layout.Column = 2;
    dds = {d1};
    if ~isempty(lab2)
        u2 = uilabel(L, 'Text', lab2); u2.Layout.Row = row; u2.Layout.Column = 3;
        d2 = uidropdown(L, 'Items', {'<none>'}, 'Value', '<none>');
        d2.Layout.Row = row; d2.Layout.Column = 4;
        dds{2} = d2;
    end
    rowOut = row + 1;
end

function setDropdown(dd, items, preferred)
    if isempty(items); items = {'<none>'}; end
    dd.Items = items;
    k = find(strcmp(items, preferred), 1);
    if ~isempty(k); dd.Value = items{k}; else; dd.Value = items{1}; end
end

function T = emptyChanTable()
    T = table((1:0)', categorical(strings(0,1), {'Neural','SlowWaveRaw','Ignore'}), ...
        strings(0,1), 'VariableNames', {'Column','Type','Label'});
end

function d = pref_dir(key)
% Starting folder for a uigetfile of type KEY: the last folder used for that
% type, else the last folder used for any file, else the current directory.
    grp = 'vengPipeline';
    d = '';
    if ispref(grp, 'lastDir'); d = getpref(grp, 'lastDir'); end
    if ispref(grp, key);       d = getpref(grp, key);       end
    if isempty(d) || ~isfolder(d); d = pwd; end
    if ~endsWith(d, filesep); d = [d filesep]; end
end

function save_pref_dir(key, folder)
% Remember FOLDER as the last location for type KEY (and as the shared last
% location). Persists across MATLAB sessions via setpref.
    if isempty(folder); return; end
    grp = 'vengPipeline';
    setpref(grp, key, folder);
    setpref(grp, 'lastDir', folder);
end

function name = guessVector(info)
    % first variable that looks like a 1-D vector
    name = info(1).name;
    for i = 1:numel(info)
        s = info(i).size;
        if numel(s) == 2 && min(s) == 1 && max(s) > 1
            name = info(i).name; return;
        end
    end
end

function name = guessMatrix(info)
    % first variable with a dimension of size 3 (the 3 antral channels),
    % else the first variable
    name = info(1).name;
    for i = 1:numel(info)
        if any(info(i).size == 3)
            name = info(i).name; return;
        end
    end
end
