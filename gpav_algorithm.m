function [YFit, vA, Blocks] = gpav_algorithm(LeftBlock, Y, W)
% GPAV_ALGORITHM  Generalised Pool Adjacent Violators for monotonic regression.
%               The GPAV algorithm follows Burdakov et al. (2006)
%   [YFit, vA, Blocks] = gpav_algorithm(LeftBlock, Y, W)
%
%   Inputs:
%       LeftBlock – cell array: for each node i, LeftBlock{i} contains indices
%                   of nodes that must be ≤ i.
%       Y         – response values (length N)
%       W         – weights (length N)
%
%   Outputs:
%       YFit  – fitted values (same order as input)
%       vA    – block‑averaged values after merging
%       Blocks – cell array of block members

    N = length(Y);
    corr = 1:N;
    vA = Y;
    xW = W;
    Blocks = cell(1,N);
    for q = 1:N
        Blocks{q} = q;
    end

    for j = 1:N
        while true
            % Find all predecessors of current block
            jMinus = LeftBlock{Blocks{j}};
            jMinus = corr(jMinus);
            Items = find(vA(jMinus) > vA(j));
            if isempty(Items)
                break
            end
            % Pick the predecessor with largest value
            LookingIn = jMinus(Items);
            [~, I] = max(vA(LookingIn));
            i = LookingIn(I(1));

            % Merge i into j
            vA(j) = (xW(i)*vA(i) + xW(j)*vA(j)) / (xW(j) + xW(i));
            xW(j) = xW(j) + xW(i);
            Blocks{j} = [Blocks{i}, Blocks{j}];
            Blocks{i} = 0;
            corr(Blocks{j}) = j;
            vA(i) = 0;
        end
    end

    YFit = vA(corr);
end