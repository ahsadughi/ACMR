function [YF, LeftBlock, Index, gtime, BlocksN] = algorithm1_basic(data, firstIncreasing, T, Trend, plotView)
% ALGORITHM1_BASIC  ACMR Algorithm 1 (basic adjacency‑constrained MR).
%
%   [YF, LeftBlock, Index, gtime, BlocksN] = algorithm1_basic(data, firstIncreasing, T, Trend, plotView)
%
%   Inputs:
%       data            – matrix with columns [X, Y, W] (time, value, weight)
%       firstIncreasing – 1 if first segment is increasing, -1 otherwise
%       T               – period length (e.g., 12 for monthly)
%       Trend           – 'D' to reverse sign for decreasing trend
%       plotView        – 0/1/2 for plotting
%
%   Outputs:
%       YF        – fitted values (original order)
%       LeftBlock – cell array of predecessors for GPAV
%       Index     – permutation vector
%       gtime     – timing vector [detect_turning, build_graph, GPAV]
%       BlocksN   – block information from GPAV

    Y = data(:,1);
    W = data(:,3);
    X = data(:,2);
    N = length(Y);
    Ksection = 2 * ceil(N / T);

    FY = Y';   % will be sign‑flipped if Trend == 'D'
    if Trend == 'D'
        FY = -FY;
        firstIncreasing = -firstIncreasing;
    end

    % Step 1: find turning points
    tic;
    [moment, ~, ~] = find_turning_points(Y, W, Ksection, firstIncreasing);
    gtime(1) = toc;

    % Step 2: build DAG (adjacency) using basic rules
    LeftBlock = cell(1,N);   % will hold indices of predecessors
    Index = zeros(1,N);
    ReIndex = zeros(1,N);
    m = 1;
    Inincreasing = firstIncreasing;
    tic;
    for j = 1:length(moment)-1
        if Inincreasing == 1
            % increasing segment
            if j == 1
                i = moment(j);
                Index(m) = i; ReIndex(i) = m; m = m+1;
            else
                i = moment(j)+1;
                if i ~= moment(j+1)+1 && i ~= moment(j+1)
                    Index(m) = i; ReIndex(i) = m;
                    if i-T > 0
                        LeftBlock{ReIndex(i)} = [ReIndex(i-T), LeftBlock{ReIndex(i)}];
                    end
                    m = m+1;
                end
            end
            i = i+1;
            upper = moment(j+1)-1;
            while i <= upper
                Index(m) = i; ReIndex(i) = m;
                m = m+1;
                LeftBlock{ReIndex(i)} = [ReIndex(i-1), LeftBlock{ReIndex(i)}];
                if i ~= upper
                    if i-T-1 > 0
                        LeftBlock{ReIndex(i)} = [ReIndex(i-T-1), LeftBlock{ReIndex(i)}];
                    end
                else
                    if i-T > 0
                        LeftBlock{ReIndex(i)} = [ReIndex(i-T), LeftBlock{ReIndex(i)}];
                    end
                end
                i = i+1;
            end
            Inincreasing = -1;
        else
            % decreasing segment
            i = moment(j+1); Lower = moment(j);
            Index(m) = i; ReIndex(i) = m;
            if i-T > 0
                LeftBlock{ReIndex(i)} = [ReIndex(i-T), LeftBlock{ReIndex(i)}];
            end
            m = m+1;
            i = i-1;
            while i >= Lower
                Index(m) = i; ReIndex(i) = m;
                LeftBlock{ReIndex(i)} = [ReIndex(i+1), LeftBlock{ReIndex(i)}];
                m = m+1;
                if i-T > 0
                    LeftBlock{ReIndex(i)} = [ReIndex(i-T), LeftBlock{ReIndex(i)}];
                end
                i = i-1;
            end
            Inincreasing = 1;
        end
    end
    % last point
    Index(N) = N; ReIndex(N) = N;
    if N-T > 0
        LeftBlock{ReIndex(N)} = [ReIndex(N-T), LeftBlock{ReIndex(N)}];
    end
    LeftBlock{ReIndex(N)} = [ReIndex(N-1), LeftBlock{ReIndex(N)}];

    gtime(2) = toc;

    % Step 3: GPAV
    vA = FY(Index);
    xW = W(Index);
    tic;
    [Yfit22, ~, BlocksN] = gpav_algorithm(LeftBlock, vA, xW);
    gtime(3) = toc;

    YF = Yfit22(Index);
    if Trend == 'D'
        YF = -YF;
    end
    YF(1) = Y(1); YF(N) = Y(N);   % fix endpoints

    % Optional plotting
    if plotView == 1
        subplot(3,1,1);
        plot(X, YF, 'linewidth',1.5); hold on;
        plot(X(moment), YF(moment), 'ro'); hold off;
        title('Fitted');
        subplot(3,1,2); plot(X, Y, 'linewidth',1.5); title('Original');
        subplot(3,1,3); plot(X, YF-Y, '*'); title('Residuals');
    elseif plotView == 2
        plot(X, Y, '-b', X, YF, '-g', 'linewidth',1.5);
        hold on; plot(X(moment), YF(moment), 'ro', 'linewidth',3); hold off;
    end
end