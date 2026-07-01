function plot_mmc(mmcFile, opts)
% PLOT_MMC  Diagnostic figures for a _mmc.mat cache written by extract_mmc.
% Three figures: (1) conditioning & cardiac QC, (2) burst detection,
% (3) metric time series. Self-contained from the cache; the raw trace in
% Fig 1A is reloaded from the source _blankmotion file if it is still findable.
%
%   plot_mmc(mmcFile)
%   plot_mmc(mmcFile, opts)
%
% opts:
%   .level   'firing' (individual muscle firings) | 'burst' (grouped episodes) ['firing']
%   .chan    channel for the time-domain / detection panels  [1]
%   .win     [t0 t1] seconds for the windowed panels          [auto: busiest 40 s]
%   .zoom    length (s) of the Fig-1C zoom inside .win         [3]
%   .saveDir folder to save the 3 figures into                [none = display only]
%   .formats subset of {'png','svg','fig'}                     [{'png','fig'}]
%
% The cache stores BOTH levels (mmc.firing and mmc.burst); .level picks which one
% the detection/metric panels display. The conditioning/cardiac QC (Fig 1) is
% level-independent. Filenames get a _firing / _burst suffix so both can be saved.

    if nargin < 2; opts = struct(); end
    if nargin < 1 || isempty(mmcFile) || exist(mmcFile,'file')~=2
        [fn,fp] = uigetfile({'*_mmc.mat','MMC cache (*_mmc.mat)'}, 'Select an _mmc.mat file');
        if isequal(fn,0); return; end
        mmcFile = fullfile(fp,fn);
    end
    S = load(mmcFile); assert(isfield(S,'mmc'), 'plot_mmc: %s has no mmc struct.', mmcFile);
    m = S.mmc; ch = getf(opts,'chan',1); zoomLen = getf(opts,'zoom',3);
    fs = m.fs; t = m.t(:); N = numel(t); sig = double(m.signal);
    col = lines(3); [~,nm] = fileparts(mmcFile);
    level = lower(getf(opts,'level','firing'));          % 'firing' | 'burst'
    lv = get_level(m, level);                            % .events/.rate/.peakAmp/.avgRate/.refractory
    Ltag = [upper(level(1)) level(2:end)];               % 'Firing' | 'Burst'

    % busiest 40 s window (most total events) unless given
    if isfield(opts,'win') && ~isempty(opts.win)
        w = opts.win;
    else
        rs = sum(lv.rate,2,'omitnan'); [~,wi] = max(rs);
        c = m.rate_t(min(max(wi,1),numel(m.rate_t)));
        w = [max(0,c-20) min(t(end),c+20)];
    end
    iw = t>=w(1) & t<=w(2);
    z0 = mean(w)-zoomLen/2; z1 = z0+zoomLen; iz = t>=z0 & t<=z1;

    % reload raw (optional)
    raw = [];
    if isfield(m,'qc') && isfield(m.qc,'srcFile') && exist(m.qc.srcFile,'file')==2
        try
            g = load(m.qc.srcFile, m.qc.dataVar); R = double(g.(m.qc.dataVar));
            if size(R,1) < size(R,2); R = R.'; end; raw = R(:, m.qc.gastricCols);
        catch; raw = []; end
    end
    rT = []; if isfield(m,'qc') && isfield(m.qc,'rpeakT'); rT = m.qc.rpeakT(:); end
    rTw = rT(rT>=w(1) & rT<=w(2));

    % ================= FIGURE 1 — conditioning & cardiac QC =================
    f1 = figure('Name',['MMC | conditioning & cardiac QC | ' nm],'Color','w','Position',[40 50 1720 980]);
    tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

    nexttile;                                            % 1A raw + R-peaks
    if ~isempty(raw)
        plot(t(iw), raw(iw,ch),'Color',[.3 .3 .3]); hold on;
        yl = ylim;
        if ~isempty(rTw)                              % all R-peak lines as ONE object
            xx = [rTw(:)'; rTw(:)'; nan(1,numel(rTw))];
            yy = repmat([yl(1); yl(2); NaN], 1, numel(rTw));
            plot(xx(:), yy(:), '-', 'Color', [1 .3 .3 .6]);
        end
        title(sprintf('1A. Raw + R-peaks (ch %d)',ch));
    else
        text(.5,.5,'raw unavailable (source file not found)','HorizontalAlignment','center'); axis off;
        title('1A. Raw + R-peaks');
    end
    xlabel('s'); xlim(w);

    nexttile;                                            % 1B conditioned + blanked gaps
    y = sig(:,ch); plot(t(iw), y(iw),'Color',col(1,:)); hold on; shade_nan(t(iw), y(iw));
    title(sprintf('1B. Conditioned (2-50 Hz) + blanked gaps (ch %d)',ch)); xlabel('s'); xlim(w);

    nexttile;                                            % 1C zoom morphology + peaks
    plot(t(iz), y(iz),'Color',col(1,:)); hold on;
    pk = find(lv.events(:,ch)); pkz = pk(t(pk)>=z0 & t(pk)<=z1);
    plot(t(pkz), y(pkz),'v','Color',col(2,:),'MarkerFaceColor',col(2,:),'MarkerSize',5);
    title(sprintf('1C. Zoom: %s morphology + peaks (%.0f-%.0f s)',level,z0,z1)); xlabel('s'); xlim([z0 z1]);

    nexttile;                                            % 1D peri-R triggered average
    if isfield(m,'qc') && ~isempty(m.qc.periR_t)
        plot(m.qc.periR_t*1e3, m.qc.periR_raw(:,ch),'Color',[.3 .3 .3]); hold on;
        plot(m.qc.periR_t*1e3, m.qc.periR_cond(:,ch),'Color',col(2,:),'LineWidth',1.5);
        legend('raw','conditioned','Location','best'); xlabel('ms re R-peak');
        title('1D. Peri-R average (conditioned should be flat)');
    else; axis off; title('1D. Peri-R average (n/a)'); end

    nexttile([1 2]);                                     % 1E PSD raw vs conditioned
    if isfield(m,'qc') && ~isempty(m.qc.psd_f)
        semilogy(m.qc.psd_f, m.qc.psd_raw(:,ch),'Color',[.3 .3 .3]); hold on;
        semilogy(m.qc.psd_f, m.qc.psd_cond(:,ch),'Color',col(2,:),'LineWidth',1.3);
        bnd = m.params.band; yl = ylim;
        plot([bnd(1) bnd(1)],yl,'k--'); plot([bnd(2) bnd(2)],yl,'k--');
        if isfield(m.qc,'meanHR') && isfinite(m.qc.meanHR)
            plot([m.qc.meanHR m.qc.meanHR]/60,yl,'r:'); % HR in Hz
        end
        legend('raw','conditioned','band','Location','best');
        xlabel('Hz'); xlim([0 100]); title('1E. PSD raw vs conditioned (slow waves removed, band kept; red = HR)');
    else; axis off; title('1E. PSD (n/a)'); end

    % ================= FIGURE 2 — burst detection ==========================
    f2 = figure('Name',['MMC | detection | ' nm],'Color','w','Position',[60 50 1720 960]);
    tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    nexttile;                                            % 2A signal + adaptive threshold + peaks
    yv = y(iw); tv = t(iw);
    win  = max(round(3*fs), round(m.params.sigmaWin*fs));
    medF = movmedian(y, win, 'omitnan');                    % full-channel moving baseline
    sgF  = movmedian(abs(y - medF), win, 'omitnan')/0.6745; % moving MAD (matches detector)
    plot(tv,yv,'Color',col(1,:)); hold on;
    plot(tv, medF(iw)+m.params.k*sgF(iw), 'r-'); plot(tv, medF(iw)-m.params.k*sgF(iw), 'r-');
    pkw = pk(t(pk)>=w(1) & t(pk)<=w(2));
    plot(t(pkw),y(pkw),'v','Color',col(2,:),'MarkerFaceColor',col(2,:),'MarkerSize',5);
    title(sprintf('2A. Conditioned + %g\\sigma threshold + %s peaks (ch %d, refr %.0f ms)', ...
        m.params.k, level, ch, lv.refractory*1e3));
    xlabel('s'); xlim(w);

    nexttile;                                            % 2B inter-event-interval hist
    pkt = t(pk); ibi = diff(pkt);
    if ~isempty(ibi); histogram(ibi, max(10,round(numel(ibi)/4))); end
    xlabel(sprintf('inter-%s interval (s)',level)); ylabel('count'); title(sprintf('2B. IEI (ch %d)',ch));

    nexttile;                                            % 2C event amplitude hist
    amp = abs(y(pk));
    if ~isempty(amp); histogram(amp, max(10,round(numel(amp)/4))); end
    xlabel(sprintf('%s peak |amplitude|',level)); ylabel('count'); title(sprintf('2C. Amplitudes (ch %d)',ch));

    nexttile; axis off;                                  % 2D QC text
    qcl = qc_lines(m, level);
    text(0.02,0.98,qcl,'VerticalAlignment','top','FontName','FixedWidth','FontSize',10,'Interpreter','none');
    title('2D. QC summary');

    % ================= FIGURE 3 — metric time series =======================
    f3 = figure('Name',['MMC | metric time series | ' nm],'Color','w','Position',[80 30 1720 1030]);
    tl = tiledlayout(5,1,'Padding','compact','TileSpacing','compact');
    ax = gobjects(5,1);

    ax(1)=nexttile;                                      % 3A rate
    plot(m.rate_t, lv.rate); hold on;
    for c3=1:3; yline(lv.avgRate(c3),'--','Color',col(c3,:)); end
    ylabel('rate (Hz)'); title(sprintf('3A. %s rate per second (3 channels; dashed = avg)',Ltag));
    legend({'ch1','ch2','ch3'},'Location','eastoutside');

    ax(2)=nexttile; plot(m.rate_t, lv.peakAmp);          % 3B peak amplitude
    ylabel('peak amp'); title(sprintf('3B. %s peak amplitude',Ltag)); legend({'ch1','ch2','ch3'},'Location','eastoutside');

    ax(3)=nexttile;                                      % 3C coverage / validity
    cov = coverage(sig, m.rate_t, m.params.W, fs)*100;
    plot(m.rate_t, cov); ylim([0 105]); ylabel('% valid');
    title('3C. Coverage per window (low = rate unreliable / NaN)'); legend({'ch1','ch2','ch3'},'Location','eastoutside');

    ax(4)=nexttile;                                      % 3D cross-channel delay
    if ~isempty(m.delay_t); plot(m.delay_t, m.delay); end
    ylabel('delay (s)'); title('3D. Inter-channel propagation delay'); legend(m.pairs,'Location','eastoutside');

    ax(5)=nexttile; hold on;                             % 3E event raster
    for c3=1:3
        et = t(lv.events(:,c3));
        plot(et, c3*ones(size(et)),'|','Color',col(c3,:),'MarkerSize',8);
    end
    ylim([0.5 3.5]); set(gca,'YTick',1:3,'YTickLabel',{'ch1','ch2','ch3'});
    xlabel('s'); title(sprintf('3E. %s event raster',Ltag));
    linkaxes(ax,'x'); xlim(ax(1),[0 t(end)]);

    % ---- match the pipeline style: white bg, black text, font 20 ----
    for ff = [f1 f2 f3]
        if exist('whiten_figure','file')==2; whiten_figure(ff); else; style_bw(ff); end
    end

    % ---- optional save ----
    saveDir = getf(opts,'saveDir','');
    if ~isempty(saveDir)
        fmts = getf(opts,'formats',{'png','fig'});
        save_figs([f1 f2 f3], {[nm '_QC_conditioning'], ...
            [nm '_QC_detection_' level], [nm '_metrics_' level]}, saveDir, fmts);
        fprintf('  [plot_mmc] saved 3 figures to %s\n', saveDir);
    end
end

% ======================================================================
function v = getf(s,f,d); if isfield(s,f)&&~isempty(s.(f)); v=s.(f); else; v=d; end; end

function lv = get_level(m, level)
% Pull the requested detection level (mmc.firing / mmc.burst). Falls back to the
% legacy single-level layout (mmc.eventSeries/rate/peakAmp/avgMMCRate) so old
% caches still plot.
    if isfield(m, level)
        lv = m.(level);
    elseif isfield(m,'eventSeries')
        rf = 0.5; if isfield(m,'params') && isfield(m.params,'refractory'); rf = m.params.refractory; end
        lv = struct('events',m.eventSeries,'rate',m.rate,'peakAmp',m.peakAmp, ...
                    'avgRate',m.avgMMCRate,'refractory',rf);
        warning('plot_mmc:legacy','legacy single-level cache; showing it as "%s". Rebuild _mmc for firing+burst.', level);
    else
        error('plot_mmc: cache has neither mmc.%s nor legacy fields -- rebuild the _mmc cache.', level);
    end
end

function style_bw(fig)
% Fallback if whiten_figure is not on the path: white bg, black text, font 20.
    FS = 20;
    set(fig,'Color','w','InvertHardcopy','off');
    for a = findall(fig,'Type','axes')'
        set(a,'Color','w','XColor','k','YColor','k','ZColor','k','GridColor',[.15 .15 .15], ...
            'FontSize',FS,'LabelFontSizeMultiplier',1,'TitleFontSizeMultiplier',1);
        try, a.Title.Color='k'; a.XLabel.Color='k'; a.YLabel.Color='k';
             a.Title.FontSize=FS; a.XLabel.FontSize=FS; a.YLabel.FontSize=FS; catch, end %#ok<CTCH>
    end
    for tlh = findall(fig,'Type','tiledlayout')'
        try, tlh.Title.Color='k'; tlh.Title.FontSize=FS; catch, end %#ok<CTCH>
    end
    for lg = findall(fig,'Type','legend')'
        try, set(lg,'Color','w','TextColor','k','EdgeColor',[.15 .15 .15],'FontSize',FS); catch, end %#ok<CTCH>
    end
    for cb = findall(fig,'Type','colorbar')'
        try, set(cb,'Color','k','FontSize',FS); cb.Label.Color='k'; cb.Label.FontSize=FS; catch, end %#ok<CTCH>
    end
    for tx = findall(fig,'Type','text')'
        try, set(tx,'FontSize',FS); c=get(tx,'Color'); if ~all(c>=0.95); set(tx,'Color','k'); end; catch, end %#ok<CTCH>
    end
end

function save_figs(figs, names, saveDir, fmts)
    if ~exist(saveDir,'dir'); mkdir(saveDir); end
    for i = 1:numel(figs)
        b = fullfile(saveDir, regexprep(names{i},'[^\w\-.]+','_'));
        if any(strcmpi(fmts,'png')); try, exportgraphics(figs(i),[b '.png'],'Resolution',150,'BackgroundColor','white'); catch, end; end %#ok<CTCH>
        if any(strcmpi(fmts,'svg')); try, exportgraphics(figs(i),[b '.svg'],'ContentType','vector','BackgroundColor','white'); catch, end; end %#ok<CTCH>
        if any(strcmpi(fmts,'fig')); try, savefig(figs(i),[b '.fig']); catch, end; end %#ok<CTCH>
    end
end

function shade_nan(tt, yy)
    bad = isnan(yy); if ~any(bad); return; end
    d = diff([0; bad(:); 0]); st = find(d==1); en = find(d==-1)-1;
    yl = ylim;
    for k=1:numel(st)
        x0=tt(st(k)); x1=tt(min(en(k),numel(tt)));
        patch([x0 x1 x1 x0],[yl(1) yl(1) yl(2) yl(2)],[1 .9 .9],'EdgeColor','none','FaceAlpha',.5);
    end
end

function cov = coverage(sig, ct, W, fs)
    N = size(sig,1); v = double(~isnan(sig)); cum = [zeros(1,3); cumsum(v,1)];
    cov = nan(numel(ct),3);
    for w=1:numel(ct)
        lo=max(1,floor((ct(w)-W/2)*fs)+1); hi=min(N,floor((ct(w)+W/2)*fs));
        cov(w,:) = (cum(hi+1,:)-cum(lo,:))/max(1,(hi-lo+1));
    end
end

function s = qc_lines(m, level)
    if nargin < 2; level = 'firing'; end
    q = []; if isfield(m,'qc'); q = m.qc; end
    L = {};
    if ~isempty(q)
        L{end+1} = sprintf('level shown : %s', level);
        L{end+1} = sprintf('R-peaks     : %d   (mean HR %.0f bpm)', numel(q.rpeakT), q.meanHR);
        L{end+1} = sprintf('%% blanked    : %s', num2str(q.pctBlanked,'%7.1f'));
        if isfield(q,'nFirings'); L{end+1} = sprintf('# firings    : %s', num2str(q.nFirings,'%7d')); end
        L{end+1} = sprintf('# bursts     : %s', num2str(q.nBursts,'%7d'));
        if isfield(m,'firing'); L{end+1} = sprintf('firings/s    : %s', num2str(m.firing.avgRate,'%7.2f')); end
        if isfield(m,'burst');  L{end+1} = sprintf('bursts/s     : %s', num2str(m.burst.avgRate, '%7.2f')); end
        if isfield(m,'firing') && isfield(m.firing,'refractory') && isfield(m,'burst')
            L{end+1} = sprintf('refractory   : %.0f ms firing / %.0f ms burst', ...
                m.firing.refractory*1e3, m.burst.refractory*1e3);
        end
        L{end+1} = sprintf('%% NaN rate   : %s', num2str(q.rateNanFrac,'%7.1f'));
        flags = {};
        if any(q.pctBlanked>40); flags{end+1}='>40% blanked'; end
        if isfield(q,'nBursts') && any(q.nBursts<10); flags{end+1}='few bursts (<10)'; end
        if any(q.rateNanFrac>50);flags{end+1}='>50% NaN rate windows'; end
        pk = abs(q.periR_cond); pr = abs(q.periR_raw);
        if any(max(pk)./max(pr) > 0.3); flags{end+1}='peri-R residual high (check cardiac removal)'; end
        if isempty(flags); flags={'none'}; end
        L{end+1} = ''; L{end+1} = ['flags: ' strjoin(flags,'; ')];
    else
        L = {'no qc fields in this cache'};
    end
    s = strjoin(L, newline);
end
