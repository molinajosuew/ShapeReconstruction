clear;
clc;

L = 10;
r = 2 ^ 7;

x = linspace(0, L, r);

[X, Y] = meshgrid(x, x);

n = 4;

for i = 0 : n
    for j = 0 : n - i
        B = nchoosek(n, i) * nchoosek(n - i, j) * X .^ i .* Y .^ j .* (L - X - Y) .^ (n - i - j) / L ^ n;
        imshow(B);
    end
end