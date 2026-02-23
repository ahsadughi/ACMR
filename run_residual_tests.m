function p = run_residual_tests(Y, Y_alg1, Y_alg2, Y_opt, Y_dp, Yraw, N, K, M, pow)
% RUN_RESIDUAL_TESTS  Apply nonparametric randomness tests to residuals.
%
%   p = run_residual_tests(Y, Y_alg1, Y_alg2, Y_opt, Y_dp, Yraw, N, K, M, pow)
%
%   Inputs: see run_simulation outputs.
%   Output: p – matrix of p‑values for each test and each scenario.

    % Load previously saved results? (Not needed if called directly)
    % For standalone use, we can loop over scenarios as in run_simulation.
    % Here we assume the inputs contain the fitted values.

    np = (K-3)*6;   % effective sample size after removing boundaries
    p = zeros(16,9);   % rows: tests, columns: scenarios

    for i = 1:9
        % Residuals from Algorithm 1
        resid1 = Yraw(2:np,i) - Y_alg1(2:np,i);
        p(1,i) = run_test(resid1);                      % run test
        [~, p(2,i)] = spearman_rho_test(resid1, 0.05); % Spearman
        [~, p(3,i)] = mann_kendall_test(resid1, 0.05); % Mann–Kendall
        [~, p(4,i)] = turning_point_test(resid1, 0.05);% turning point

        % Residuals from Algorithm 2
        resid2 = Yraw(2:np,i) - Y_alg2(2:np,i);
        p(5,i) = run_test(resid2);
        [~, p(6,i)] = spearman_rho_test(resid2, 0.05);
        [~, p(7,i)] = mann_kendall_test(resid2, 0.05);
        [~, p(8,i)] = turning_point_test(resid2, 0.05);

        % Residuals from optimal solution
        resid_opt = Yraw(2:np,i) - Y_opt(2:np,i);
        p(9,i) = run_test(resid_opt);
        [~, p(10,i)] = spearman_rho_test(resid_opt, 0.05);
        [~, p(11,i)] = mann_kendall_test(resid_opt, 0.05);
        [~, p(12,i)] = turning_point_test(resid_opt, 0.05);

        % Residuals from DP
        resid_dp = Yraw(2:np,i) - Y_dp(2:np,i);
        p(13,i) = run_test(resid_dp);
        [~, p(14,i)] = spearman_rho_test(resid_dp, 0.05);
        [~, p(15,i)] = mann_kendall_test(resid_dp, 0.05);
        [~, p(16,i)] = turning_point_test(resid_dp, 0.05);
    end
end