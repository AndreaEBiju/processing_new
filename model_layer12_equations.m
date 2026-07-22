function varargout = model_layer12_equations(what, varargin)
%MODEL_LAYER12_EQUATIONS Dispatcher for Layer 1-2 model building blocks.
%   Usage: model_layer12_equations('h', u_E, u_half)
%          model_layer12_equations('Phi_c', u_E, f_max_c)
%          model_layer12_equations('g1', u_M, f_lo, f_hi)
%          model_layer12_equations('g2', u_M, u_M_half)
%          model_layer12_equations('f1', u_E, u_half, f_max_c)
%          model_layer12_equations('f2', u_E, u_half, f_max_c2)
%          model_layer12_equations('r_vagus', u_M, u_E, params)  % params = struct
%
% See model_fitting_plan.md Section 2 for equations/justification and Section 4
% for known limitations (L1-L11). Do not add terms or "improve" these equations
% without flagging back to the user first (per that document's Section 6).
%
% params (for 'r_vagus') is a struct with fields:
%   W_c, u_half, f_max_c, f_lo, f_hi, u_M_half, f_max_c2, w_1, w_2
% These field names are canonical for this task -- reuse exactly, do not rename.
%
% Two corrections from an earlier draft, both load-bearing (see
% model_fitting_plan.md Sections 2.1-2.3 and matlab_implementation_instructions.md
% Section 3):
%   1. w_c and lambda_max_cap are symbolically non-separable (they multiply
%      together everywhere either appears) -- merged into a single parameter W_c.
%      There is no longer a separate 'lambda_cap' dispatcher case; it was folded
%      into term_c directly inside r_vagus_fn.
%   2. f2(u_E) = Phi_c(u_E) alone does NOT vanish at u_E=0 (Phi_c(0)=1), which
%      silently violates the "no standalone mechanical term" requirement (L1).
%      Fixed by gating f2 with h(u_E) (like f1) but giving it its OWN following
%      corner f_max_c2, distinct from f1's f_max_c, so f1 and f2 stay
%      non-proportional (required for the moving-optimum property, Section 3).

    switch what
        % 'term_shares' is the only multi-output case (returns [p_cap,p_1,p_2]);
        % everything else returns exactly one value into varargout{1}.
        case 'term_shares'
            [varargout{1:max(nargout,1)}] = term_shares_fn(varargin{:});
            return;
        case 'h';        out = h_fn(varargin{:});
        case 'Phi_c';     out = Phi_c_fn(varargin{:});
        case 'g1';        out = g1_fn(varargin{:});
        case 'g2';        out = g2_fn(varargin{:});
        case 'f1';        out = f1_fn(varargin{:});
        case 'f2';        out = f2_fn(varargin{:});
        case 'r_vagus';   out = r_vagus_fn(varargin{:});
        case 'cv2_hat';   out = cv2_hat_fn(varargin{:});
        case 'fwhm_hat';  out = fwhm_hat_fn(varargin{:});
        otherwise; error('model_layer12_equations:unknown','Unknown request: %s', what);
    end
    varargout{1} = out;
end

function y = h_fn(u_E, u_half)
% Electrical entrainment, saturating, no threshold (model_fitting_plan.md Sec 2.1).
    y = u_E ./ (u_E + u_half);
end

function y = Phi_c_fn(u_E, f_max_c)
% C-fibre following. f_max_c is the CALLER-SUPPLIED corner -- pass p.f_max_c for
% term_c/term_1's corner, or p.f_max_c2 for term_2's own, distinct corner. Do not
% hardcode which corner this uses; it is a generic function of whichever corner
% is passed. NOTE: Phi_c(0) = 1, not 0 -- this is correct for Phi_c in isolation.
% It is f2 (below), gated by h(u_E), that must vanish at u_E=0, not Phi_c itself.
    y = f_max_c ./ (f_max_c + u_E);
end

function y = g1_fn(u_M, f_lo, f_hi)
% Bandpass mechanical coupling, tuned toward the M50xE100 hotspot (Sec 2.2).
% L4: this shape is NOT independently validated at receptor level -- justified
% only by the interaction pattern, not controlled single-modality dose-response
% data. L10: paired with g2 below, this is a rank-2 mathematical convenience for
% ONE unresolved mechanoreceptor population, not two distinct receptor types.
    y = (u_M ./ (u_M + f_lo)) .* (f_hi ./ (u_M + f_hi));
end

function y = g2_fn(u_M, u_M_half)
% Monotonic saturating mechanical coupling, tuned toward the M100xE10 hotspot
% (Sec 2.2). See L4/L10 note on g1_fn above -- applies here too.
    y = u_M ./ (u_M + u_M_half);
end

function y = f1_fn(u_E, u_half, f_max_c)
    y = h_fn(u_E, u_half) .* Phi_c_fn(u_E, f_max_c);
end

function y = f2_fn(u_E, u_half, f_max_c2)
% CORRECTED: previously Phi_c_fn(u_E, f_max_c) alone, which does NOT vanish at
% u_E=0 (Phi_c(0)=1) -- violated the no-standalone-mechanical-term requirement
% (L1). Now gated by h(u_E) like f1, but with its OWN corner f_max_c2 (must
% differ from f1's f_max_c, else f1 and f2 become proportional and the model
% loses the ability to shift its preferred u_M with u_E -- see
% model_fitting_plan.md Section 3).
    y = h_fn(u_E, u_half) .* Phi_c_fn(u_E, f_max_c2);
end

function y = r_vagus_fn(u_M, u_E, p)
% p is a struct with fields: W_c, u_half, f_max_c, f_lo, f_hi, u_M_half,
%                            f_max_c2, w_1, w_2
% NOTE: W_c replaces the old separate w_c/lambda_max_cap (merged -- see above).
    term_c  = p.W_c .* Phi_c_fn(u_E, p.f_max_c) .* h_fn(u_E, p.u_half);
    term_1  = p.w_1 .* f1_fn(u_E, p.u_half, p.f_max_c)  .* g1_fn(u_M, p.f_lo, p.f_hi);
    term_2  = p.w_2 .* f2_fn(u_E, p.u_half, p.f_max_c2) .* g2_fn(u_M, p.u_M_half);
    y = term_c + term_1 + term_2;
    % L1: NO standalone mechanical term -- every term above now correctly
    % vanishes at u_E=0 by construction (term_c/term_1 via the shared h(u_E)
    % factor, term_2 via its own h(u_E) factor). This is a parsimony choice
    % under identifiability constraints (9 params vs 9 M x E conditions, L9),
    % NOT a claim that the true standalone mechanical effect is proven zero.
    % Verified explicitly in test_layer12_equations.m -- this property was
    % violated in an earlier draft and must not regress silently.
end

function [p_cap, p_1, p_2] = term_shares_fn(u_M, u_E, p)
% model_fitting_plan.md Section 2.4: term-share mixture weights, shared by
% cv2_hat/fwhm_hat. Same three underlying shapes (T_cap/T_1/T_2) as r_vagus_fn,
% renormalized to sum to 1 -- but T_1/T_2 here use v_1/v_2, NOT rate's w_1/w_2.
%
% ARCHITECTURE CHANGE (per user, after real-data finding that w_1/w_2 fit to
% ~0 for rate): reusing rate's w_1/w_2 here forces p_1/p_2 ~= 0 everywhere
% whenever the rate interaction is null, making it IMPOSSIBLE to represent
% "rate synergy ~ 0, but CV2/FWHM synergy real" -- confirmed on real data
% (p_cap=1.000000, p_1/p_2 ~1e-18 at every M x E condition, CV2/FWHM Stage 3
% fit showed 0.00% improvement over a constant-mean baseline as a direct
% mechanical consequence, not evidence of no CV2/FWHM effect). v_1/v_2 are a
% SEPARATE free parameter pair, specific to the CV2/FWHM mixture, fit in
% fit_layer12_stage3_shape.m -- deliberately named v_1/v_2 (not w_1_cv2 etc.)
% so it is visually unambiguous in code/logs that these are distinct from
% rate's w_1/w_2, not a variant of them. L2 still applies: this whole mixture
% approach is a first-pass placeholder, not derived from spike-train statistics.
    T_cap = p.W_c .* Phi_c_fn(u_E, p.f_max_c) .* h_fn(u_E, p.u_half);
    T_1   = p.v_1 .* f1_fn(u_E, p.u_half, p.f_max_c)  .* g1_fn(u_M, p.f_lo, p.f_hi);
    T_2   = p.v_2 .* f2_fn(u_E, p.u_half, p.f_max_c2) .* g2_fn(u_M, p.u_M_half);
    total = T_cap + T_1 + T_2;
    % Guard the u_M=0 & u_E=0 corner where all three terms are exactly zero by
    % construction (L1) -- shares are undefined there (0/0), not zero; NaN is
    % the honest answer, not silently reported as e.g. all-cap. This condition
    % is never one of the fitted M x E combined trials (those all have
    % u_M>0 & u_E>0), so it does not affect Stage 3 fitting, only direct calls.
    safe = total ~= 0;
    p_cap = nan(size(total)); p_1 = p_cap; p_2 = p_cap;
    p_cap(safe) = T_cap(safe) ./ total(safe);
    p_1(safe)   = T_1(safe)   ./ total(safe);
    p_2(safe)   = T_2(safe)   ./ total(safe);
end

function y = cv2_hat_fn(u_M, u_E, p)
% model_fitting_plan.md Section 2.4 / Section 5.4. p extends the r_vagus params
% with CV2_cap, CV2_1, CV2_2 (the "pure" CV2 if each term totally dominated).
% Linear-mixture approximation, NOT derived from actual ISI statistics (L2).
    [p_cap, p_1, p_2] = term_shares_fn(u_M, u_E, p);
    y = p_cap.*p.CV2_cap + p_1.*p.CV2_1 + p_2.*p.CV2_2;
end

function y = fwhm_hat_fn(u_M, u_E, p)
% model_fitting_plan.md Section 2.4 / Section 5.4. p extends the r_vagus params
% with FWHM_cap, FWHM_1, FWHM_2. See cv2_hat_fn above -- identical structure.
    [p_cap, p_1, p_2] = term_shares_fn(u_M, u_E, p);
    y = p_cap.*p.FWHM_cap + p_1.*p.FWHM_1 + p_2.*p.FWHM_2;
end
