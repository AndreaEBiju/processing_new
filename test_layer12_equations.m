classdef test_layer12_equations < matlab.unittest.TestCase
%TEST_LAYER12_EQUATIONS Unit tests for model_layer12_equations.m.
%   Run via: runtests('test_layer12_equations')
%   run_layer12_first_pass.m aborts before fitting if any of these fail --
%   do not fit against real data with failing unit tests (see
%   matlab_implementation_instructions.md Section 7 checklist).

    methods (Test)
        function test_r_vagus_zero_at_uE_zero(tc)
            % REGRESSION TEST: an earlier draft of f2 used Phi_c(u_E) alone,
            % which does NOT vanish at u_E=0 (Phi_c(0)=1) and silently broke
            % this exact property. This test must pass on the CURRENT f2
            % (gated by h(u_E) with its own corner f_max_c2) -- if it starts
            % failing again, f2 has regressed to the old broken form.
            p = default_test_params();
            uM = [10 50 100];
            for i = 1:numel(uM)
                y = model_layer12_equations('r_vagus', uM(i), 0, p);
                tc.verifyEqual(y, 0, 'AbsTol', 1e-9);
            end
        end

        function test_f2_zero_at_uE_zero_specifically(tc)
            % Isolates f2 itself (not just the full r_vagus sum) to make the
            % regression this guards against unambiguous if it ever reappears.
            y = model_layer12_equations('f2', 0, 50, 100);
            tc.verifyEqual(y, 0, 'AbsTol', 1e-9);
        end

        function test_f1_f2_not_proportional(tc)
            % f_max_c and f_max_c2 must differ, else f1 and f2 collapse to
            % rescaled copies of each other and the model loses the
            % moving-optimum property (model_fitting_plan.md Section 3).
            % Checks that the ratio f1/f2 is NOT constant across u_E.
            uE = [10 50 100 500 1000];
            u_half = 50; f_max_c = 200; f_max_c2 = 60;  % deliberately distinct
            f1v = model_layer12_equations('f1', uE, u_half, f_max_c);
            f2v = model_layer12_equations('f2', uE, u_half, f_max_c2);
            ratio = f1v ./ f2v;
            tc.verifyTrue(range(ratio) > 1e-6, ...
                'f1/f2 ratio is constant -- f_max_c and f_max_c2 must differ.');
        end

        function test_h_monotonic_increasing(tc)
            uE = linspace(0.1, 2000, 500);
            y = model_layer12_equations('h', uE, 100);
            tc.verifyTrue(all(diff(y) >= 0));
            tc.verifyTrue(all(y >= 0 & y < 1));
        end

        function test_Phi_c_monotonic_decreasing(tc)
            uE = linspace(0.1, 2000, 500);
            y = model_layer12_equations('Phi_c', uE, 100);
            tc.verifyTrue(all(diff(y) <= 0));
            tc.verifyTrue(all(y > 0 & y <= 1));
            % NOTE: Phi_c(0) = 1, not 0 -- this is correct and expected for
            % Phi_c in isolation. It is f2 (gated by h) that must vanish at
            % u_E=0, not Phi_c itself. Do not "fix" this test to expect 0.
        end

        function test_g1_has_single_interior_peak(tc)
            uM = linspace(0.1, 500, 1000);
            y = model_layer12_equations('g1', uM, 20, 80);
            [~, idx] = max(y);
            tc.verifyTrue(idx > 1 && idx < numel(uM));      % interior, not monotonic
            tc.verifyTrue(all(diff(y(1:idx)) >= -1e-9));    % rising before peak
            tc.verifyTrue(all(diff(y(idx:end)) <= 1e-9));   % falling after peak
        end

        function test_g2_monotonic_increasing_bounded(tc)
            uM = linspace(0.1, 500, 500);
            y = model_layer12_equations('g2', uM, 50);
            tc.verifyTrue(all(diff(y) >= 0));
            tc.verifyTrue(all(y >= 0 & y < 1));
        end

        function test_no_division_errors_at_zero(tc)
            p = default_test_params();
            tc.verifyWarningFree(@() model_layer12_equations('r_vagus', 0, 0, p));
            tc.verifyWarningFree(@() model_layer12_equations('r_vagus', 0, 100, p));
        end

        function test_r_vagus_positive_for_combined_conditions(tc)
            % Sanity check: with all-positive params and u_M,u_E > 0, the
            % compound signal should be strictly positive (all three terms
            % are products of positive saturating/bandpass functions).
            p = default_test_params();
            [uMg, uEg] = meshgrid([10 50 100], [10 100 1000]);
            y = model_layer12_equations('r_vagus', uMg, uEg, p);
            tc.verifyTrue(all(y(:) > 0));
        end

        function test_r_vagus_vectorized_matches_elementwise(tc)
            % Vectorized call must match looped scalar calls exactly -- guards
            % against accidental non-elementwise operators creeping in.
            p = default_test_params();
            uM = [10 50 100 10 50 100 10 50 100];
            uE = [10 10 10 100 100 100 1000 1000 1000];
            yVec = model_layer12_equations('r_vagus', uM, uE, p);
            yLoop = arrayfun(@(m,e) model_layer12_equations('r_vagus', m, e, p), uM, uE);
            tc.verifyEqual(yVec, yLoop, 'AbsTol', 1e-12);
        end

        function test_term_shares_sum_to_one(tc)
            % model_fitting_plan.md Section 2.4: p_cap+p_1+p_2 must sum to 1
            % everywhere the three terms aren't all exactly zero.
            p = default_test_params_shape();
            [uMg, uEg] = meshgrid([10 50 100], [10 100 1000]);
            [p_cap, p_1, p_2] = model_layer12_equations('term_shares', uMg, uEg, p);
            tc.verifyEqual(p_cap + p_1 + p_2, ones(size(uMg)), 'AbsTol', 1e-10);
        end

        function test_term_shares_nan_at_uM0_uE0(tc)
            % At u_M=0 & u_E=0, T_cap=T_1=T_2=0 by construction (L1) -- shares
            % are undefined (0/0), must be NaN, not silently 0 or all-cap.
            p = default_test_params_shape();
            [p_cap, p_1, p_2] = model_layer12_equations('term_shares', 0, 0, p);
            tc.verifyTrue(isnan(p_cap) && isnan(p_1) && isnan(p_2));
        end

        function test_term_shares_use_v1_v2_not_w1_w2(tc)
            % ARCHITECTURE CHECK: term_shares must use v_1/v_2, NOT w_1/w_2 --
            % changing w_1/w_2 alone (v_1/v_2 held fixed) must NOT change the
            % shares. This guards against silently reverting to the old
            % (rate-coupled) wiring that made p_1/p_2 track rate's w_1/w_2.
            p = default_test_params_shape();
            uM = 50; uE = 100;
            [p_cap1, p_1_1, p_2_1] = model_layer12_equations('term_shares', uM, uE, p);
            p2 = p; p2.w_1 = p.w_1 * 37; p2.w_2 = p.w_2 * 0.001;   % change rate-only weights
            [p_cap2, p_1_2, p_2_2] = model_layer12_equations('term_shares', uM, uE, p2);
            tc.verifyEqual([p_cap1,p_1_1,p_2_1], [p_cap2,p_1_2,p_2_2], 'AbsTol', 1e-12);
            % but changing v_1/v_2 (rate-independent) MUST change the shares
            p3 = p; p3.v_1 = p.v_1 * 5; p3.v_2 = p.v_2 * 5;
            [~, p_1_3, ~] = model_layer12_equations('term_shares', uM, uE, p3);
            tc.verifyNotEqual(p_1_3, p_1_1);
        end

        function test_cv2_fwhm_hat_bounded_by_reference_values(tc)
            % CV2_hat/FWHM_hat are convex combinations (weights sum to 1, all
            % nonnegative since every term_share component is a product of
            % nonnegative functions) of the three reference values -- so they
            % must lie within [min(refs), max(refs)] everywhere shares are defined.
            p = default_test_params_shape();
            [uMg, uEg] = meshgrid([10 50 100], [10 100 1000]);
            cv2v = model_layer12_equations('cv2_hat', uMg, uEg, p);
            fwhmv = model_layer12_equations('fwhm_hat', uMg, uEg, p);
            cv2Refs = [p.CV2_cap, p.CV2_1, p.CV2_2];
            fwhmRefs = [p.FWHM_cap, p.FWHM_1, p.FWHM_2];
            tc.verifyTrue(all(cv2v(:) >= min(cv2Refs) - 1e-9 & cv2v(:) <= max(cv2Refs) + 1e-9));
            tc.verifyTrue(all(fwhmv(:) >= min(fwhmRefs) - 1e-9 & fwhmv(:) <= max(fwhmRefs) + 1e-9));
        end

        function test_cv2_hat_matches_manual_computation(tc)
            % Concrete numeric check at one (u_M,u_E) point against a
            % hand-computed term-share mixture, independent of term_shares_fn.
            p = default_test_params_shape();
            uM = 50; uE = 100;
            T_cap = p.W_c * model_layer12_equations('Phi_c', uE, p.f_max_c) * model_layer12_equations('h', uE, p.u_half);
            T_1   = p.v_1 * model_layer12_equations('f1', uE, p.u_half, p.f_max_c) * model_layer12_equations('g1', uM, p.f_lo, p.f_hi);
            T_2   = p.v_2 * model_layer12_equations('f2', uE, p.u_half, p.f_max_c2) * model_layer12_equations('g2', uM, p.u_M_half);
            total = T_cap + T_1 + T_2;
            expected = (T_cap/total)*p.CV2_cap + (T_1/total)*p.CV2_1 + (T_2/total)*p.CV2_2;
            y = model_layer12_equations('cv2_hat', uM, uE, p);
            tc.verifyEqual(y, expected, 'AbsTol', 1e-10);
        end
    end
end

function p = default_test_params()
    p = struct('W_c',1,'u_half',50,'f_max_c',200, ...
                'f_lo',20,'f_hi',80,'u_M_half',50,'f_max_c2',60, ...
                'w_1',1,'w_2',1);
    % f_max_c2 deliberately != f_max_c (200 vs 60) -- see
    % test_f1_f2_not_proportional above for why this must hold.
end

function p = default_test_params_shape()
    p = default_test_params();
    % v_1/v_2 deliberately != w_1/w_2 -- these are the SEPARATE, CV2/FWHM-only
    % mixture weights (see term_shares_fn header); tests must exercise that
    % term_shares/cv2_hat/fwhm_hat use v_1/v_2, never w_1/w_2.
    p.v_1 = 2; p.v_2 = 3;
    p.CV2_cap = 0.3; p.CV2_1 = 0.9; p.CV2_2 = 0.5;      % deliberately distinct
    p.FWHM_cap = 0.8; p.FWHM_1 = 1.6; p.FWHM_2 = 1.1;   % deliberately distinct
end
