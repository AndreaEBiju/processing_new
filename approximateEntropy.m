function ApEn = approximateEntropy(rr, m, r)
% approximateEntropy  Compute Approximate Entropy (ApEn)
%
%   ApEn = approximateEntropy(rr, m, r)
%
%   rr  - time series
%   m   - embedding dimension (typically 2)
%   r   - tolerance (e.g., 0.15–0.25 * std(rr))

    if nargin < 2, m = 2; end
    rr = rr(:)';
    N  = length(rr);
    if nargin < 3, r = 0.2 * std(rr); end

    phi = zeros(1,2);

    for mm = m:(m+1)
        Cm = zeros(1, N-mm+1);

        for i = 1:(N-mm+1)
            xmi = rr(i:(i+mm-1));
            count = 0;

            for j = 1:(N-mm+1)
                xmj = rr(j:(j+mm-1));
                if max(abs(xmi - xmj)) <= r
                    count = count + 1;
                end
            end

            Cm(i) = count / (N-mm+1);
        end

        phi(mm-m+1) = mean(log(Cm));
    end

    ApEn = phi(1) - phi(2);
end
