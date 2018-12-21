function NSP = NonSepBern(L, r, s, n, i, j)
    [x_1, x_2] = meshgrid(linspace(0, L, r), linspace(0, L, s));
    NSP = L ^ (- n) * nchoosek(n, i) * nchoosek(n - i , j) * x_1 .^ i .* x_2 .^ j .* (L - x_1 - x_2) .^ (n - i - j);
end