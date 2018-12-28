function NSBM = NonSepBernMon(L, a, b, n, i, j)
    [X, Y] = meshgrid(linspace(0, L, a), linspace(0, L, b));
    NSBM = zeros(a, b);
    for p = 0 : n - i - j
        for q = 0 : n - i - j - p
            NSBM = NSBM + nchoosek(n, i) * nchoosek(n - i, j) * factorial(n - i - j) / (factorial(p) * factorial(q) * factorial(n - i - j - p - q)) * L ^ (p - n) * (- 1) ^ (n - i - j - p) * X .^ (q + i) .* Y .^ (n - i - j - p - q + j);
        end
    end
end