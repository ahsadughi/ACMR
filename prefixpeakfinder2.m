function Verror = prefixpeakfinder2(y, w)
% PREFIXPEAKFINDER2  Cumulative forward isotonic errors for a unimodal sequence.
%
%   Verror = prefixpeakfinder2(y, w)
%
%   Inputs:
%       y – vector of data values
%       w – vector of weights (same length)
%
%   Output:
%       Verror – cumulative squared errors after isotonic regression
%                from the beginning up to each index.
%
%   Used by find_turning_points.m to locate peaks.

    y = [-Inf; y(:)];
    w = [0; w(:)];
    n = length(y);
    Vmean = zeros(n,1);
    Vleft = zeros(n,1);
    sumwy = zeros(n,1);
    sumwy2 = zeros(n,1);
    sumw = zeros(n,1);
    Verror = zeros(n,1);

    Vmean(1) = -Inf;
    Vleft(1) = 0;
    Verror(1) = 0;

    for i = 2:n
        Vmean(i) = y(i);
        Vleft(i) = i;
        sumwy(i) = w(i) * y(i);
        sumwy2(i) = w(i) * y(i)^2;
        sumw(i) = w(i);

        % Merge while the current block mean is ≤ the previous block mean
        while (i > 1) && (Vmean(i) <= Vmean(Vleft(i)-1))
            % Merge with previous block
            sumwy(i) = w(i)*y(i) + w(Vleft(i)-1)*y(Vleft(i)-1);
            sumwy2(i) = w(i)*y(i)^2 + w(Vleft(i)-1)*y(Vleft(i)-1)^2;
            sumw(i) = w(i) + w(Vleft(i)-1);
            Vmean(i) = sumwy(i) / sumw(i);
            Vleft(i) = Vleft(Vleft(i)-1);
        end

        levelerror = sumwy2(i) - sumwy(i)^2 / sumw(i);
        if i ~= 1
            Verror(i) = levelerror + Verror(i-1);
        end
    end

    % Remove the artificial first element
    Verror = Verror(2:n);
end