function NSBM = NonSepBernMonomial(L, a, b, n, i, j)
[Y, X] = meshgrid(linspace(0, L, b), linspace(0, L, a));
NSBM = zeros(a, b);
for p = 0 : n - i - j
    for q = 0 : n - i - j - p
        NSBM = NSBM + nchoosek(n, i) * nchoosek(n - i, j) * factorial(n - i - j) / (factorial(p) * factorial(q) * factorial(n - i - j - p - q)) * L ^ (p - n) * (- 1) ^ (n - i - j - p) * X .^ (q + i) .* Y .^ (n - i - j - p - q + j);
    end
end
end