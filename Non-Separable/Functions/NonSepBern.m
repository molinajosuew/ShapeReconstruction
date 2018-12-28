function NSB = NonSepBern(L, a, b, n, i, j)
    [x, y] = meshgrid(linspace(0, L, a), linspace(0, L, b));
    NSB = L ^ (- n) * nchoosek(n, i) * nchoosek(n - i , j) * x .^ i .* y .^ j .* (L - x - y) .^ (n - i - j);
end