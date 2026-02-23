function [YF, adj, Index, gtime, BlocksN] = algorithm2_transfer(data, firstIncreasing, T, Trend, plotView)
% ALGORITHM2_TRANSFER  ACMR Algorithm 2 (transfer‑function version).
%
%   [YF, adj, Index, gtime, BlocksN] = algorithm2_transfer(data, firstIncreasing, T, Trend, plotView)
%
%   Inputs: see algorithm1_basic.
%   Outputs:
%       adj     – sparse adjacency matrix
%       Index   – permutation vector (reordering used)
%       BlocksN – block information from GPAV

    Y = data(:,1);
    W = data(:,3);
    X = data(:,2);
    N = length(Y);
    Ksection = 2 * ceil(N / T);

    FY = Y';
    if Trend == 'D'
        FY = -FY;
        firstIncreasing = -firstIncreasing;
    end

    % Step 1: find turning points
    tic;
    [moment, ~, ~] = find_turning_points(Y, W, Ksection, firstIncreasing);
    gtime(1) = toc;

    % Step 2: build adjacency matrix with transfer function
    adj = sparse(N, N);
    Index = zeros(1,N);
    ReIndex = zeros(1,N);
    m = 1;
    Inincreasing = firstIncreasing;
    tic;
    for j = 1:length(moment)-1
        mm = X(moment(j));
        Mx = X(moment(j+1));
        d = Mx - mm;
        if Inincreasing == 1
            % increasing segment
            if j == 1
                i = moment(j);
                Index(m) = i; ReIndex(i) = m; m = m+1;
            else
                i = moment(j)+1;
                if i ~= moment(j+1)+1 && i ~= moment(j+1)
                    Index(m) = i; ReIndex(i) = m;
                    % transfer‑function edge
                    if j > 2
                        u = (X(i)-mm)/d;
                        l = u*(X(moment(j-1))-X(moment(j-2))) + X(moment(j-2));
                        hb = find(X(moment(j-2):moment(j-1)) <= l, 1, 'last');
                        hb = hb + moment(j-2) - 1;
                        adj(ReIndex(hb), ReIndex(i)) = 1;
                    end
                    m = m+1;
                end
            end
            i = i+1;
            upper = moment(j+1)-1;
            while i <= upper
                Index(m) = i; ReIndex(i) = m; m = m+1;
                adj(ReIndex(i-1), ReIndex(i)) = 1;
                if j > 2
                    u = (X(i)-mm)/d;
                    l = u*(X(moment(j-1))-X(moment(j-2))) + X(moment(j-2));
                    hb = find(X(moment(j-2):moment(j-1)) <= l, 1, 'last');
                    hb = hb + moment(j-2) - 1;
                    adj(ReIndex(hb), ReIndex(i)) = 1;
                end
                i = i+1;
            end
            Inincreasing = -1;
        else
            % decreasing segment
            i = moment(j+1); Lower = moment(j);
            Index(m) = i; ReIndex(i) = m;
            if j > 2
                u = (X(i)-mm)/d;
                l = u*(X(moment(j-1))-X(moment(j-2))) + X(moment(j-2));
                hb = find(X(moment(j-2):moment(j-1)) >= l, 1);
                hb = hb + moment(j-2) - 1;
                adj(ReIndex(hb), ReIndex(i)) = 1;
            end
            m = m+1;
            i = i-1;
            while i >= Lower
                Index(m) = i; ReIndex(i) = m;
                adj(ReIndex(i+1), ReIndex(i)) = 1;
                if j > 2
                    u = (X(i)-mm)/d;
                    l = u*(X(moment(j-1))-X(moment(j-2))) + X(moment(j-2));
                    hb = find(X(moment(j-2):moment(j-1)) >= l, 1);
                    hb = hb + moment(j-2) - 1;
                    adj(ReIndex(hb), ReIndex(i)) = 1;
                end
                m = m+1;
                i = i-1;
            end
            Inincreasing = 1;
            if moment(j)-1 > 0
                adj(ReIndex(moment(j)-1), ReIndex(moment(j))) = 1;
            end
        end
    end
    % last point
    Index(N) = N; ReIndex(N) = N;
    if N-T > 0
        adj(ReIndex(N-T), ReIndex(N)) = 1;
    end

    gtime(2) = toc;

    % Step 3: GPAV
    vA = FY(Index);
    xW = W(Index);
    tic;
    [Yfit22, ~, BlocksN] = gpav_algorithm(adj, vA, xW);  % note: uses adjacency matrix
    gtime(3) = toc;

    YF = Yfit22(Index);
    if Trend == 'D'
        YF = -YF;
    end
    YF(1) = Y(1); YF(N) = Y(N);

    % Optional plotting (same as algorithm1_basic)
    % ... (omitted for brevity)
end