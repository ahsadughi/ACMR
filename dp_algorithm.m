function [YFIT, G, ITAU, ITHETA, sumnb] = dp_algorithm(data, Ksection, first, plotFlag)
% DP_ALGORITHM  Dynamic programming solution for piecewise monotonic regression.
%
%   [YFIT, G, ITAU, ITHETA, sumnb] = dp_algorithm(data, Ksection, first, plotFlag)
%
%   Inputs:
%       data     – matrix with columns [value, time, weight] (order may differ)
%       Ksection – number of sections (turning points)
%       first    – 'I' for increasing first segment, 'D' for decreasing
%       plotFlag – 0/1 to plot results
%
%   Outputs:
%       YFIT   – fitted values
%       G      – accumulated cost matrix
%       ITAU   – index of optimal path
%       ITHETA – estimated turning points
%       sumnb  – number of blocks (degrees of freedom)

    F = data(:,1);          % values
    WF = data(:,3);         % weights
    X = data(:,2);          % time

    % Find local extrema
    [LOWER, IUPPER, Index, NL, NU] = funcXTREMA(F, Ksection);
    N = length(F);
    I1 = 1;

    % Determine turning points ITHETA
    ITHETA = zeros(1, Ksection);
    IL = I1; IU = I1;
    for J = 1:Ksection
        if mod(J,2) ~= 0
            ITHETA(J) = LOWER(IL);
            IL = IL + 1;
        else
            ITHETA(J) = IUPPER(IU);
            IU = IU + 1;
        end
    end
    IL = I1; IU = I1;

    % If decreasing trend, flip sign
    if first == 'D'
        F = -F;
    end

    % Prepare reversed arrays
    FT = zeros(1,N);
    WFT = zeros(1,N);
    FTNEG = zeros(1,N);
    for I = I1:N
        FT(I) = F(N + I1 - I);
        WFT(I) = WF(N + I1 - I);
        FTNEG(I) = -F(N + I1 - I);
    end

    % Initialise G and ITAU
    G = zeros(Ksection, Ksection);
    ITAU = zeros(Ksection, Ksection);
    [SS, ~, ~] = funcL2WMON(F, WF, first, I1, N);
    for J = I1:NU
        G(1, J) = SS(IUPPER(J));
        ITAU(1, J) = I1;
    end

    % Main dynamic programming loops
    ILU = IUPPER(IU);
    while ILU <= LOWER(NL-1)
        IL = IL + 1;
        ILU = LOWER(IL);
        [SS, ~, ~] = funcL2WMON(FT, WFT, first, N+I1-ILU, N+I1-ITHETA(2));
        for M = 2:2:Ksection-1
            if M == 2
                lim = I1;
            else
                lim = ITAU(M-2, IL);
            end
            if lim < ITHETA(M)
                K = Index(ITHETA(M));
            else
                K = Index(lim);
            end
            RG = zeros(1, Ksection);
            for J = K:IU
                RG(J) = G(M-1, J) + SS(N + I1 - IUPPER(J));
            end
            [~, MINRG] = min(RG(K:IU));
            MINRG = MINRG + K - 1;
            G(M, IL) = RG(MINRG);
            ITAU(M, IL) = IUPPER(MINRG);
            ITHETA(M) = IUPPER(MINRG);
        end

        if ILU > IUPPER(NU-1)
            break
        end
        IU = IU + 1;
        ILU = IUPPER(IU);
        [SS, ~, ~] = funcL2WMON(FTNEG, WFT, first, N+I1-ILU, N+I1-ITHETA(3));
        for M = 3:2:Ksection-1
            lim = ITAU(M-2, IU);
            if lim < ITHETA(M)
                K = Index(ITHETA(M));
            else
                K = Index(lim);
            end
            RG = zeros(1, Ksection);
            for J = K:IL
                RG(J) = G(M-1, J) + SS(N + I1 - LOWER(J));
            end
            [~, MINRG] = min(RG(K:IL));
            MINRG = MINRG + K - 1;
            G(M, IU) = RG(MINRG);
            ITAU(M, IU) = LOWER(MINRG);
            ITHETA(M) = LOWER(MINRG);
        end
    end

    % Final segment
    if mod(Ksection,2) == 0
        if Ksection > 2
            K = ITAU(Ksection-2, IL);
        else
            K = 2;
        end
        [SS, ~, ~] = funcL2WMON(FT, WFT, first, I1, N+I1-K);
        K = Index(K);
        RG = zeros(1, NU);
        for J = K:NU
            RG(J) = G(Ksection-1, J) + SS(N + I1 - IUPPER(J));
        end
        [~, MINRG] = min(RG(K:NU));
        MINRG = MINRG + K - 1;
        G(Ksection, NU) = RG(MINRG);
        ITAU(Ksection, NU) = IUPPER(MINRG);
        ITHETA(Ksection) = IUPPER(MINRG);
    else
        K = ITAU(Ksection-2, IU);
        [SS, ~, ~] = funcL2WMON(FT, WFT, first, I1, N+I1-K);
        K = Index(K);
        RG = zeros(1, NL-1);
        for J = K:NL-1
            RG(J) = G(Ksection-1, J) + SS(N + I1 - LOWER(J));
        end
        [~, MINRG] = min(RG(K:NL-1));
        MINRG = MINRG + K - 1;
        G(Ksection, NL) = RG(MINRG);
        ITAU(Ksection, NL) = LOWER(MINRG);
        ITHETA(Ksection) = LOWER(MINRG);
    end

    % Backtrack to recover all turning points
    if mod(Ksection,2) == 0
        MINRG = NU;
        ITHETA(Ksection) = N;
        if ITHETA(Ksection-1) == 0
            MINRG = Index(I1);
        elseif Ksection > 2
            MINRG = Index(ITHETA(Ksection-1));
        else
            MINRG = Index(I1+1);
        end
        for J = Ksection-2:-1:1
            ITHETA(J) = ITAU(J+1, MINRG);
            if ITHETA(J) == 0
                MINRG = Index(I1);
            else
                MINRG = Index(ITHETA(J));
            end
        end
        if Ksection > 2
            ITHETA(Ksection-2) = ITHETA(Ksection-1);
        else
            I1 = ITHETA(Ksection-1);
        end
        ITHETA(Ksection-1) = ITHETA(Ksection);
    else
        MINRG = NL;
        ITHETA(Ksection) = N;
        if ITHETA(Ksection-1) == 0
            MINRG = Index(I1);
        elseif Ksection > 3
            MINRG = Index(ITHETA(Ksection-3));
        else
            MINRG = Index(I1);
        end
        for J = Ksection-4:-1:1
            ITHETA(J) = ITAU(J+1, MINRG);
            if ITHETA(J) == 0
                MINRG = Index(I1);
            else
                MINRG = Index(ITHETA(J));
            end
        end
        if Ksection > 3
            ITHETA(Ksection-3) = ITHETA(Ksection-2);
        else
            ITHETA(I1) = ITHETA(Ksection-2);
        end
        ITHETA(Ksection-2) = ITHETA(Ksection-1);
        ITHETA(Ksection-1) = ITHETA(Ksection);
    end

    % Reconstruct fitted values
    sumnb = 0;
    YFIT = zeros(1, N);
    for i = 1:Ksection-1
        if i == 1
            L1 = I1;
        else
            L1 = ITHETA(i-1) + 1;
        end
        LN = ITHETA(i);
        if i == Ksection-2
            LN = N;
        end
        if mod(i,2) == 0
            [SS, Yseg, ~, nb] = funcL2WMON(FT, WFT, 1, L1, LN);
            sumnb = sumnb + nb;
            YFIT(L1:LN) = Yseg;
        else
            [SS, FTNEGseg, ~, nb] = funcL2WMON(FT, WFT, 1, N+I1-LN, N+I1-L1);
            sumnb = sumnb + nb;
            for I = L1:LN
                YFIT(I) = FTNEG(N + I1 - I);
            end
        end
    end
    % Ensure exact fit at turning points
    for i = 1:Ksection
        YFIT(ITHETA(i)) = F(ITHETA(i));
    end

    if plotFlag
        figure;
        subplot(2,1,1); plot(X, YFIT); title('Fitted');
        subplot(2,1,2); plot(X, F); title('Original');
    end
end