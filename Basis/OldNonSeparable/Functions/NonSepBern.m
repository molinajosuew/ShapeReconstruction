function NSB = NonSepBern(L, a, b, n, i, j)
[X, Y] = meshgrid(linspace(0, L, b), linspace(0, L, a));
NSB = L ^ (- n) * nchoosek(n, i) * nchoosek(n - i , j) * X .^ i .* Y .^ j .* (L - X - Y) .^ (n - i - j);
end