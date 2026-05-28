%% Demo: boxScatterPlot with synthetic data
clc; clear; close all;

rng(42);

% Pick a theme — 'light' or 'dark'. boxScatterPlot inherits the palette
% and background that pubfig_setup installs.
pubfig_setup('Theme','light');
% pubfig_setup('Theme','dark');    % uncomment to preview dark mode

outDir = fullfile(pwd,'plots_demo');
if ~exist(outDir,'dir'), mkdir(outDir); end

%% (1) Ungrouped: 3 conditions, N=20 per condition
d1 = { 72 + 4*randn(20,1), ...     % cond1: baseline HR
       68 + 5*randn(20,1), ...     % cond2: stim
       70 + 4.5*randn(20,1) };     % cond3: recovery

f1 = figure('Position',[100 100 700 480],'Name','demo_ungrouped');
boxScatterPlot(d1, ...
    'GroupLabels', {'baseline','stim','recovery'}, ...
    'YLabel',      'HR (bpm)', ...
    'Title',       'Heart rate across conditions');
exportgraphics(f1, fullfile(outDir,'demo_ungrouped.png'), 'Resolution',150);

%% (2) Grouped: 3 conditions x 3 time windows, paired across windows
nSub = 8;
condShift = [0 -3 -1];           % per-condition mean shift
winShift  = [0 -1.5 -3];         % per-window mean shift (longer window -> more drop)

d2 = cell(3,3);
for g = 1:3
    subj = 72 + 4*randn(nSub,1);     % subject-specific baseline
    for s = 1:3
        d2{g,s} = subj + condShift(g) + winShift(s) + 0.8*randn(nSub,1);
    end
end

f2 = figure('Position',[100 100 900 520],'Name','demo_grouped');
boxScatterPlot(d2, ...
    'GroupLabels',    {'cond1','cond2','cond3'}, ...
    'SubgroupLabels', {'1 min','2 min','5 min'}, ...
    'YLabel',         'HR (bpm)', ...
    'Title',          'Heart rate by condition and window', ...
    'ColorBySubject', true);
exportgraphics(f2, fullfile(outDir,'demo_grouped.png'), 'Resolution',150);

%% (3) Grouped without pairing (independent samples per window)
d3 = cell(3,3);
for g = 1:3
    for s = 1:3
        d3{g,s} = 72 + condShift(g) + winShift(s) + 4*randn(18,1);
    end
end

f3 = figure('Position',[100 100 900 520],'Name','demo_grouped_unpaired');
boxScatterPlot(d3, ...
    'GroupLabels',    {'cond1','cond2','cond3'}, ...
    'SubgroupLabels', {'1 min','2 min','5 min'}, ...
    'YLabel',         'HR (bpm)', ...
    'Title',          'Grouped, unpaired samples');
exportgraphics(f3, fullfile(outDir,'demo_grouped_unpaired.png'), 'Resolution',150);

fprintf('Saved demos to: %s\n', outDir);
