% This script corroborates that the expansion coefficients work.

clear;
clc;

n = 4;
k = 4;
L = 500;

N = ProtoNonSeparableBernsteinExpansionCoefficients(n, L);
N = N(:, :, :, 501 : end, 501 : end);
[X, Y] = meshgrid(0 : L, 0 : L);

for i = 0 : k
    for j = 0 : k - i
        I1 = conv2(wjmBSpline(k, 1), wjmBSpline(k, 1), reshape(N(1 + k, 1 + i, 1 + j, :, :), size(N, 4), size(N, 5)), "valid");
        I2 = nchoosek(k, i) .* nchoosek(k - i, j) .* X .^ i .* Y .^ j .* (L - X - Y) .^ (k - i - j) .* L .^ (- k);
        
        surf(abs(I1 - I2));
        
        pause(.5);
    end
end