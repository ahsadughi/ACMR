function [moment, INDEX, Kstep] = find_turning_points(Y, W, Ksection, firstIncreasing)
% FIND_TURNING_POINTS  Locate turning points via prefix isotonic regression.
%
%   [moment, INDEX, Kstep] = find_turning_points(Y, W, Ksection, firstIncreasing)
%
%   Inputs:
%       Y, W           – data and weights
%       Ksection       – desired number of sections (turning points)
%       firstIncreasing – 1 if first segment increasing, -1 otherwise
%
%   Outputs:
%       moment – vector of turning point indices
%       INDEX  – (unused) placeholder
%       Kstep  – approximate step size

    N = length(Y);
    Kstep = ceil(N / Ksection);
    j = 1;
    lomax = firstIncreasing;
    i = 1;
    moment = [];

    while i < N
        if i == 1
            Kstepone = Kstep - 1;
            sf = prefixpeakfinder2(lomax * Y(i:i+Kstepone), W(i:i+Kstepone));
            sb = prefixpeakfinder2(lomax * Y(i+Kstepone:-1:i), W(i+Kstepone:-1:i));
            sb = sb(Kstepone+1:-1:1);
        elseif i <= N - Kstep
            sf = prefixpeakfinder2(lomax * Y(i:i+Kstep-1), W(i:i+Kstep-1));
            sb = prefixpeakfinder2(lomax * Y(i+Kstep-1:-1:i), W(i+Kstep-1:-1:i));
            sb = sb(Kstep:-1:1);
        else
            sf = prefixpeakfinder2(lomax * Y(i:N), W(i:N));
            sb = prefixpeakfinder2(lomax * Y(N:-1:i), W(N:-1:i));
            sb = sb(N-i+1:-1:1);
        end
        lomax = -lomax;
        [~, idx] = min(sf + sb);
        moment(j) = idx + i - 1;
        j = j + 1;
        i = i + Kstep;
    end

    if moment(1) ~= 1
        moment = [1, moment];
    end
    moment = [moment, N];
end

% Note: prefixpeakfinder2 is not provided; it should compute the isotonic
% regression and return the cumulative squared errors.