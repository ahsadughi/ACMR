function [x, resnorm, cpuTime] = optimal_solution(Index, adj, data, MaxIter)
% OPTIMAL_SOLUTION  Solve the constrained quadratic program.
%
%   [x, resnorm, cpuTime] = optimal_solution(Index, adj, data, MaxIter)
%
%   Inputs:
%       Index   – permutation vector (from algorithm2_transfer)
%       adj     – adjacency matrix (sparse logical)
%       data    – raw data matrix [X, Y, W]
%       MaxIter – maximum iterations for lsqlin
%
%   Outputs:
%       x       – optimal fitted values (in original order)
%       resnorm – squared norm of residuals
%       cpuTime – elapsed CPU time

    % Reorder adjacency according to Index
    adj2 = adj(Index, Index);
    N = length(adj);
    [I, J] = find(adj2);
    nCon = nnz(I);
    Aineq = zeros(nCon, N);
    for k = 1:nCon
        Aineq(k, I(k)) = -1;
        Aineq(k, J(k)) = 1;
    end
    bineq = zeros(nCon, 1);

    C = eye(N);
    d = data(:,1);   % original Y

    options = optimset('MaxIter', MaxIter, 'Display', 'off', 'LargeScale', 'off');
    tic;
    [x, resnorm] = lsqlin(C, d, Aineq, bineq, [], [], [], [], [], options);
    cpuTime = toc;
end