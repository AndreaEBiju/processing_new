function [alpha1, alpha2, nVals, F] = dfaRR(rr, nVals)
% dfaRR  Detrended Fluctuation Analysis for RR series
%
%   [alpha1, alpha2, nVals, F] = dfaRR(rr, nVals)
%
%   rr     - time series (RR intervals)
%   nVals  - (optional) vector of box sizes (in samples)
%
%   alpha1 - short-term scaling exponent (e.g. 4–16 beats)
%   alpha2 - long-term scaling exponent (>16 beats)
%   nVals  - box sizes used
%   F      - fluctuation function values

    x = rr(:) - mean(rr);
    N = length(x);

    % Integrated profile
    y = cumsum(x);

    % Default box sizes
    if nargin < 2
        nMin = 4;
        nMax = floor(N/4);
        nVals = unique(round(logspace(log10(nMin), log10(nMax), 15)));
    end

    F = zeros(length(nVals),1);

    for k = 1:length(nVals)
        n = nVals(k);
        Ns = floor(N / n);
        if Ns < 2
            F(k) = NaN;
            continue;
        end

        % Truncate to multiple of n
        yTrunc = y(1:Ns*n);
        yMat = reshape(yTrunc, n, Ns);

        % Detrend each segment
        t = (1:n)';
        F2 = zeros(Ns,1);
        for v = 1:Ns
            p = polyfit(t, yMat(:,v), 1);
            trend = polyval(p, t);
            F2(v) = mean((yMat(:,v) - trend).^2);
        end

        F(k) = sqrt(mean(F2));
    end

    % Remove NaNs
    valid = ~isnan(F) & F>0;
    nVals = nVals(valid);
    F     = F(valid);

    logn = log10(nVals);
    logF = log10(F);

    % Short-term region: 4–16 beats
    shortIdx = nVals >= 4 & nVals <= 16;
    if sum(shortIdx) >= 3
        p1 = polyfit(logn(shortIdx), logF(shortIdx), 1);
        alpha1 = p1(1);
    else
        alpha1 = NaN;
    end

    % Long-term region: > 16 beats
    longIdx = nVals > 16;
    if sum(longIdx) >= 3
        p2 = polyfit(logn(longIdx), logF(longIdx), 1);
        alpha2 = p2(1);
    else
        alpha2 = NaN;
    end
end
