function [results, Yraw] = run_simulation(N, K, M, pow)
% RUN_SIMULATION  Main simulation driver for ACMR paper.
%
%   [results, Yraw] = run_simulation(N, K, M, pow)
%
%   Inputs:
%       N   - sample size
%       K   - number of periods (roughly N/12 for monthly data)
%       M   - number of replications for timing measurements
%       pow - trend type: 1 for linear, ~0.7 for nonlinear
%
%   Outputs:
%       results - matrix with objective values, CPU times, relative errors
%       Yraw    - raw simulated data for each error distribution (9 columns)
%
%   The function generates data for 9 scenarios (3 error distributions × 3 noise levels),
%   applies Algorithm 1, Algorithm 2, DP, and the optimal quadratic program,
%   and collects performance metrics.

    T = 12;   % period length (months)

    % Generate data for all nine scenarios
    [rawdataNL, rawdataNM, rawdataNH, ...
     rawdataEL, rawdataEM, rawdataEH, ...
     rawdataDL, rawdataDM, rawdataDH] = funcSimulaDataGen6(N, pow);

    Yraw = [rawdataNL(:,1), rawdataNM(:,1), rawdataNH(:,1), ...
            rawdataEL(:,1), rawdataEM(:,1), rawdataEH(:,1), ...
            rawdataDL(:,1), rawdataDM(:,1), rawdataDH(:,1)];

    nn = length(rawdataNM);

    % Preallocate output arrays
    T_alg1 = zeros(1,9);   R_alg1 = zeros(1,9);   R2_alg1 = zeros(1,9);
    T_alg2 = zeros(1,9);   R_alg2 = zeros(1,9);   R2_alg2 = zeros(1,9);
    T_dp   = zeros(1,9);   R_dp   = zeros(1,9);   R2_dp   = zeros(1,9);
    T_opt  = zeros(1,9);   R_opt  = zeros(1,9);   R2_opt  = zeros(1,9);

    Y_alg1 = zeros(nn,9);
    Y_alg2 = zeros(nn,9);
    Y_dp   = zeros(nn,9);
    Y_opt  = zeros(nn,9);

    % Loop over the nine scenarios
    for s = 1:9
        switch s
            case 1, data = rawdataNL;
            case 2, data = rawdataNM;
            case 3, data = rawdataNH;
            case 4, data = rawdataEL;
            case 5, data = rawdataEM;
            case 6, data = rawdataEH;
            case 7, data = rawdataDL;
            case 8, data = rawdataDM;
            case 9, data = rawdataDH;
        end

        % ---------- Algorithm 1 (basic) ----------
        tic;
        for rep = 1:M
            [YF, ~, ~, ~, BlocksN] = algorithm1_basic(data, 1, T, 'D', 0);
        end
        T_alg1(s) = toc / M;
        R_alg1(s) = sum((data(:,1) - YF').^2);
        Y_alg1(:,s) = YF';
        nb = count_blocks(BlocksN, length(YF));
        R2_alg1(s) = sum((data(:,1) - YF').^2) / (length(YF) - 1.5*nb);

        % ---------- Algorithm 2 (transfer) ----------
        tic;
        for rep = 1:M
            [YF, adj, Index, ~, ~] = algorithm2_transfer(data, 1, T, 'D', 0);
        end
        T_alg2(s) = toc / M;
        R_alg2(s) = sum((data(:,1) - YF').^2);
        Y_alg2(:,s) = YF';
        % (Blocks not available from algorithm2_transfer, skip adj MSE)

        % ---------- DP algorithm ----------
        tic;
        for rep = 1:M
            [YFIT, ~, ~, ~, sumnb] = dp_algorithm(data, K, 'I', 0);
        end
        T_dp(s) = toc / M;
        R_dp(s) = sum((data(1:nn-24,1) - YFIT(1:nn-24)').^2);  % DP discards last 24?
        Y_dp(:,s) = YFIT;
        nb = min(length(YFIT)/3, sumnb);
        R2_dp(s) = sum((data(:,1) - YFIT').^2) / (length(YFIT) - 1.5*nb);

        % ---------- Optimal solution (quadratic programming) ----------
        % First run algorithm2_transfer to get adj and Index needed for exact
        [~, adj, Index] = algorithm2_transfer(data, 1, T, 'D', 0);
        tic;
        [x, resnorm, ~] = optimal_solution(Index, adj, data, 100000);
        T_opt(s) = toc;
        R_opt(s) = sum((data(2:nn-1,1) - x(2:nn-1)).^2);   % interior points only
        Y_opt(2:nn-1,s) = x(2:nn-1);
        Y_opt(1,s) = data(1,1);
        Y_opt(nn,s) = data(nn,1);
        R2_opt(s) = resnorm / (length(adj) - 1);
    end

    % Collect results in a matrix (as in original paper)
    results = zeros(9,13);
    results(:,1)  = R_alg1';
    results(:,2)  = (R_alg1 - R_opt) ./ R_opt;
    results(:,3)  = (R_alg1 - R_dp) ./ R_dp;
    results(:,4)  = T_alg1';
    results(:,5)  = R_opt';
    results(:,6)  = (R_opt - R_dp) ./ R_dp;
    results(:,7)  = T_opt';
    results(:,8)  = R_dp';
    results(:,9)  = T_dp';
    results(:,10) = N;
    results(:,11) = K;
    results(:,12) = M;
    results(:,13) = pow;

end

% Helper: count non‑empty blocks
function nb = count_blocks(BlocksN, len)
    nb = 0;
    for i = 1:len
        if ~isempty(BlocksN{i}) && BlocksN{i} ~= 0
            nb = nb + 1;
        end
    end
    nb = min(len/3, nb);
end