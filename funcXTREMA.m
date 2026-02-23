function [LOWER, IUPPER, INDX, NL, NU] = funcXTREMA(F, KSECTN)
% FUNCXTREMA  Find indices of local minima (LOWER) and maxima (IUPPER).
%
%   [LOWER, IUPPER, INDX, NL, NU] = funcXTREMA(F, KSECTN)
%
%   Inputs:
%       F      – vector of data values (length N)
%       KSECTN – desired number of sections (turning points)
%
%   Outputs:
%       LOWER  – indices of local minima (strictly decreasing then increasing)
%       IUPPER – indices of local maxima (strictly increasing then decreasing)
%       INDX   – mapping from original index to its position in LOWER/IUPPER
%       NL     – number of elements in LOWER
%       NU     – number of elements in IUPPER
%
%   Implements the algorithm described in Demetriou & Powell (1991).

    I1 = 1;
    N = size(F,1);
    LOWER = zeros(1,N);   % pre‑allocate
    IUPPER = zeros(1,N);
    INDX = zeros(1,N);

    IL = I1;
    LOWER(IL) = I1;
    INDX(I1) = IL;
    IU = I1 - 1;

    % Find first plateau
    IA = I1;
    while (IA < N) && (F(IA) == F(IA+1))
        IA = IA + 1;
    end

    % Trivial case: all equal
    if IA == N
        IU = I1;
        IUPPER(IU) = I1;
        INDX(I1) = IU;
        if mod(KSECTN,2) == 0
            IL = IL + 1;
            LOWER(IL) = N;
            INDX(N) = IL;
        else
            IU = IU + 1;
            IUPPER(IU) = N;
            INDX(N) = IU;
        end
        NL = IL;
        NU = IU;
        return;
    end

    % Include I1 in IUPPER if it is a peak start
    if F(IA) > F(IA+1)
        IU = I1;
        IUPPER(IU) = I1;
        INDX(I1) = IU;
    end

    % Find last plateau
    IB = N;
    while (IB > 1) && (F(IB) == F(IB-1))
        IB = IB - 1;
    end

    % Scan interior points
    for I = IA+1 : IB-1
        II = I;
        while (II < IB-1) && (F(II) == F(II+1))
            II = II + 1;
        end
        % Local minimum
        if (F(I-1) > F(I)) && (F(II) < F(II+1))
            IL = IL + 1;
            LOWER(IL) = I;
            INDX(I) = IL;
        end
        % Local maximum
        if (F(I-1) < F(I)) && (F(II) > F(II+1))
            IU = IU + 1;
            IUPPER(IU) = I;
            INDX(I) = IU;
        end
    end

    % Include N in LOWER or IUPPER according to parity and slope
    if (mod(KSECTN,2) == 0) || (F(IB-1) > F(IB))
        IL = IL + 1;
        LOWER(IL) = N;
        INDX(N) = IL;
    end
    if (mod(KSECTN,2) ~= 0) || (F(IB-1) < F(IB))
        IU = IU + 1;
        IUPPER(IU) = N;
        INDX(N) = IU;
    end

    NL = IL;
    NU = IU;
end