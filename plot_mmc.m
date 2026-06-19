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
%   .chan   channel for the time-domain / detection panels   [1]
%   .win    [t0 t1] seconds for the windowed panels          [auto: busiest 40 s]
%   .zoom   length (s) of the Fig-1C zoom inside .win         [3]

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

    % busiest 40 s window (most total bursts) unless given
    if isfield(opts,'win') && ~isempty(opts.win)
        w = opts.win;
    else
        rs = sum(m.rate,2,'omitnan'); [~,wi] = max(rs);
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
    figure('Name',['MMC | conditioning & cardiac QC | ' nm],'Color','w','Position',[60 60 1180 760]);
    tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

    nexttile;                                            % 1A raw + R-peaks
    if ~isempty(raw)
        plot(t(iw), raw(iw,ch),'Color',[.3 .3 .3]); hold on;
        yl = ylim; for k=1:numel(rTw); plot([rTw(k) rTw(k)],yl,'r-','Color',[1 .3 .3 .5]); end
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
    pk = find(m.eventSeries(:,ch)); pkz = pk(t(pk)>=z0 & t(pk)<=z1);
    plot(t(pkz), y(pkz),'v','Color',col(2,:),'MarkerFaceColor',col(2,:),'MarkerSize',5);
    title(sprintf('1C. Zoom: burst morphology + peaks (%.0f-%.0f s)',z0,z1)); xlabel('s'); xlim([z0 z1]);

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
    figure('Name',['MMC | detection | ' nm],'Color','w','Position',[90 90 1180 720]);
    tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    nexttile;                                            % 2A signal + adaptive threshold + peaks
    yv = y(iw); tv = t(iw);
    win = max(3,round(m.params.sigmaWin*fs));
    med = movmedian(yv,win,'omitnan'); sg = movmedian(abs(yv-med),win,'omitnan')/0.6745;
    plot(tv,yv,'Color',col(1,:)); hold on;
    plot(tv,med+m.params.k*sg,'r-'); plot(tv,med-m.params.k*sg,'r-');
    pkw = pk(t(pk)>=w(1) & t(pk)<=w(2));
    plot(t(pkw),y(pkw),'v','Color',col(2,:),'MarkerFaceColor',col(2,:),'MarkerSize',5);
    title(sprintf('2A. Conditioned + %g\\sigma threshold + burst peaks (ch %d)',m.params.k,ch));
    xlabel('s'); xlim(w);

    nexttile;                                            % 2B inter-burst-interval hist
    pkt = t(pk); ibi = diff(pkt);
    if ~isempty(ibi); histogram(ibi, max(10,round(numel(ibi)/4))); end
    xlabel('inter-burst interval (s)'); ylabel('count'); title(sprintf('2B. IBI (ch %d)',ch));

    nexttile;                                            % 2C burst amplitude hist
    amp = abs(y(pk));
    if ~isempty(amp); histogram(amp, max(10,round(numel(amp)/4))); end
    xlabel('burst peak |amplitude|'); ylabel('count'); title(sprintf('2C. Amplitudes (ch %d)',ch));

    nexttile; axis off;                                  % 2D QC text
    qcl = qc_lines(m);
    text(0.02,0.98,qcl,'VerticalAlignment','top','FontName','FixedWidth','FontSize',10,'Interpreter','none');
    title('2D. QC summary');

    % ================= FIGURE 3 — metric time series =======================
    figure('Name',['MMC | metric time series | ' nm],'Color','w','Position',[120 60 1180 820]);
    tl = tiledlayout(5,1,'Padding','compact','TileSpacing','compact');
    ax = gobjects(5,1);

    ax(1)=nexttile;                                      % 3A rate
    plot(m.rate_t, m.rate); hold on;
    for c3=1:3; yline(m.avgMMCRate(c3),'--','Color',col(c3,:)); end
    ylabel('rate (Hz)'); title('3A. Burst firing rate per second (3 channels; dashed = avg)');
    legend({'ch1','ch2','ch3'},'Location','eastoutside');

    ax(2)=nexttile; plot(m.rate_t, m.peakAmp);           % 3B peak amplitude
    ylabel('peak amp'); title('3B. Burst peak amplitude'); legend({'ch1','ch2','ch3'},'Location','eastoutside');

    ax(3)=nexttile;                                      % 3C coverage / validity
    cov = coverage(sig, m.rate_t, m.params.W, fs)*100;
    plot(m.rate_t, cov); ylim([0 105]); ylabel('% valid');
    title('3C. Coverage per window (low = rate unreliable / NaN)'); legend({'ch1','ch2','ch3'},'Location','eastoutside');

    ax(4)=nexttile;                                      % 3D cross-channel delay
    if ~isempty(m.delay_t); plot(m.delay_t, m.delay); end
    ylabel('delay (s)'); title('3D. Inter-channel propagation delay'); legend(m.pairs,'Location','eastoutside');

    ax(5)=nexttile; hold on;                             % 3E event raster
    for c3=1:3
        et = t(m.eventSeries(:,c3));
        plot(et, c3*ones(size(et)),'|','Color',col(c3,:),'MarkerSize',8);
    end
    ylim([0.5 3.5]); set(gca,'YTick',1:3,'YTickLabel',{'ch1','ch2','ch3'});
    xlabel('s'); title('3E. Burst event raster');
    linkaxes(ax,'x'); xlim(ax(1),[0 t(end)]);
end

% ======================================================================
function v = getf(s,f,d); if isfield(s,f)&&~isempty(s.(f)); v=s.(f); else; v=d; end; end

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

function s = qc_lines(m)
    q = []; if isfield(m,'qc'); q = m.qc; end
    L = {};
    if ~isempty(q)
        L{end+1} = sprintf('R-peaks   : %d   (mean HR %.0f bpm)', numel(q.rpeakT), q.meanHR);
        L{end+1} = sprintf('%% blanked  : %s', num2str(q.pctBlanked,'%6.1f'));
        L{end+1} = sprintf('# bursts   : %s', num2str(q.nBursts,'%6d'));
        L{end+1} = sprintf('avg rate Hz: %s', num2str(m.avgMMCRate,'%6.3f'));
        L{end+1} = sprintf('%% NaN rate : %s', num2str(q.rateNanFrac,'%6.1f'));
        flags = {};
        if any(q.pctBlanked>40); flags{end+1}='>40% blanked'; end
        if any(q.nBursts<10);    flags{end+1}='few bursts (<10)'; end
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
