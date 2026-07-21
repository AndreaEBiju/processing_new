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
    end
end

function p = default_test_params()
    p = struct('W_c',1,'u_half',50,'f_max_c',200, ...
                'f_lo',20,'f_hi',80,'u_M_half',50,'f_max_c2',60, ...
                'w_1',1,'w_2',1);
    % f_max_c2 deliberately != f_max_c (200 vs 60) -- see
    % test_f1_f2_not_proportional above for why this must hold.
end
