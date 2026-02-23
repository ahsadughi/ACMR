function [rawdataNL, rawdataNM, rawdataNH, ...
          rawdataEL, rawdataEM, rawdataEH, ...
          rawdataDL, rawdataDM, rawdataDH] = funcSimulaDataGen6(N, pow)
% FUNCSIMULADATAGEN6  Generate simulated data for the ACMR paper.
%
%   [rawdataNL, rawdataNM, rawdataNH, ...] = funcSimulaDataGen6(N, pow)
%
%   Inputs:
%       N   – desired sample size (approximate, actual length may differ)
%       pow – exponent for the trend component: 1 for linear, 0.7 for nonlinear
%
%   Outputs (each is a matrix with columns [value, time, weight]):
%       rawdataNL – Normal error, low  noise
%       rawdataNM – Normal error, medium noise
%       rawdataNH – Normal error, high  noise
%       rawdataEL – Exponential error, low  noise
%       rawdataEM – Exponential error, medium noise
%       rawdataEH – Exponential error, high  noise
%       rawdataDL – Double‑exponential (Laplace) error, low  noise
%       rawdataDM – Double‑exponential error, medium noise
%       rawdataDH – Double‑exponential error, high  noise
%
%   The underlying model is  M(t) = sin(2π t + t^0.005) - t^pow.
%   The first 12 observations are replaced by the noise‑free signal to
%   avoid edge effects (as in the paper).

    NN = ceil(N/12);
    t = 0:(1/12):NN+0;
    x = sin(2*pi*t + t.^0.005);          % seasonal component with slight phase shift
    dim1 = length(t);
    dim2 = dim1 - 8;                      % effective length after trimming

    % --- Normal errors ---
    R = normrnd(0, 0.5, 1, length(t));    % low noise
    NL = x - t.^pow - R;
    R = normrnd(0, 1.0, 1, length(t));    % medium noise
    NM = x - t.^pow - R;
    R = normrnd(0, 1.2, 1, length(t));    % high noise
    NH = x - t.^pow - R;

    rawdataNL = [(NL(9:dim1))', (1:dim2)', ones(dim2,1)];
    rawdataNM = [(NM(9:dim1))', (1:dim2)', ones(dim2,1)];
    rawdataNH = [(NH(9:dim1))', (1:dim2)', ones(dim2,1)];

    % --- Exponential errors ---
    s = 0.5;   r = -s * log(rand(1,length(t)));   EL = x - t.^pow - r;
    s = 1.0;   r = -s * log(rand(1,length(t)));   EM = x - t.^pow - r;
    s = 1.2;   r = -s * log(rand(1,length(t)));   EH = x - t.^pow - r;

    rawdataEL = [(EL(9:dim1))', (1:dim2)', ones(dim2,1)];
    rawdataEM = [(EM(9:dim1))', (1:dim2)', ones(dim2,1)];
    rawdataEH = [(EH(9:dim1))', (1:dim2)', ones(dim2,1)];

    % --- Double‑exponential (Laplace) errors ---
    z = rand(length(t),1);
    zz = zeros(length(t),1);
    in = z <= 0.5;   ip = z > 0.5;

    lambda = 5.0;    % low noise
    zz(in) =  1/lambda * log(2*z(in));
    zz(ip) = -1/lambda * log(2*(1-z(ip)));
    DL = x - t.^pow - zz';

    lambda = 2.5;    % medium noise
    zz(in) =  1/lambda * log(2*z(in));
    zz(ip) = -1/lambda * log(2*(1-z(ip)));
    DM = x - t.^pow - zz';

    lambda = 1.0;    % high noise
    zz(in) =  1/lambda * log(2*z(in));
    zz(ip) = -1/lambda * log(2*(1-z(ip)));
    DH = x - t.^pow - zz';

    rawdataDL = [(DL(9:dim1))', (1:dim2)', ones(dim2,1)];
    rawdataDM = [(DM(9:dim1))', (1:dim2)', ones(dim2,1)];
    rawdataDH = [(DH(9:dim1))', (1:dim2)', ones(dim2,1)];

    % Replace first 12 observations with the noise‑free signal (as in paper)
    tt = x - t.^pow;
    rawdataNL(1:12,1) = tt(1:12)';
    rawdataNM(1:12,1) = tt(1:12)';
    rawdataNH(1:12,1) = tt(1:12)';
    rawdataEL(1:12,1) = tt(1:12)';
    rawdataEM(1:12,1) = tt(1:12)';
    rawdataEH(1:12,1) = tt(1:12)';
    rawdataDL(1:12,1) = tt(1:12)';
    rawdataDM(1:12,1) = tt(1:12)';
    rawdataDH(1:12,1) = tt(1:12)';
end