clear;
clc;

n = 4;
L = 500;
N = ProtoNonSeparableBernsteinExpansionCoefficients(n, L);

[X, Y] = meshgrid(0 : L, 0 : L);

for i = 0 : n
    for j = 0 : n - i
        I1 = conv2(wjmBSpline(n, 1), wjmBSpline(n, 1), reshape(N(1 + i, 1 + j, :, :), size(N, 3), size(N, 4)), "valid");
        I2 = nchoosek(n, i) .* nchoosek(n - i, j) .* X .^ i .* Y .^ j .* (L - X - Y) .^ (n - i - j) .* L .^ (- n);
        
        imshow(abs(I1 - I2));
        
        pause(.5);
    end
end