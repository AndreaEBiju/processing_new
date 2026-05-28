function sampEn = sampleEntropy(rr, m, r)
% sampleEntropy  Compute Sample Entropy (SampEn) of RR series
%
%   sampEn = sampleEntropy(rr, m, r)
%
%   rr  - time series (e.g. RR intervals)
%   m   - embedding dimension (typically 2)
%   r   - tolerance (e.g. 0.15–0.25 * std(rr))
%
%   sampEn - Sample entropy value

    if nargin < 2, m = 2; end
    rr = rr(:)';
    N  = length(rr);
    if nargin < 3, r = 0.2 * std(rr); end

    % Count matches for m and m+1
    B = 0; A = 0;  % B: for m, A: for m+1

    for i = 1:(N - m)
        xmi = rr(i:(i+m-1));
        xm1i = rr(i:(i+m));

        for j = (i+1):(N - m)
            xmj  = rr(j:(j+m-1));
            xm1j = rr(j:(j+m));

            if max(abs(xmi - xmj)) < r
                B = B + 1;
                if max(abs(xm1i - xm1j)) < r
                    A = A + 1;
                end
            end
        end
    end

    if B == 0 || A == 0
        sampEn = NaN;   % no matches -> undefined
    else
        sampEn = -log(A / B);
    end
end
