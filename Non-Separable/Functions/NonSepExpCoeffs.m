function NSEC = NonSepExpCoeffs(n, a, b)
C1 = DualExpansionCoefficients(n, a, 1);
C2 = DualExpansionCoefficients(n, b, 1);
NSEC = zeros(n + 1, n + 1, size(C2, 2), size(C1, 2));
for i = 0 : n
    for j = 0 : n - i
        for q = 0 : n - i - j
            for r = 0 : n - i - j - q
                NSEC(1 + i, 1 + j, :, :) = NSEC(1 + i, 1 + j, :, :) + reshape(C1(1 + r + j, :)' * nchoosek(n, i) * nchoosek(n - i, j) * factorial(n - i - j) / (factorial(n - i - j - q - r) * factorial(q) * factorial(r)) * (- 1) ^ (q + r) * C2(1 + q + i, :), 1, 1, size(C2, 2), size(C1, 2));
            end
        end
    end
end
end