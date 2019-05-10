function NSEC = NonSepExpCoeffs(n, r, s)
C1 = DualExpansionCoefficients(n, r, 1);
C2 = DualExpansionCoefficients(n, s, 1);
NSEC = zeros(n + 1, n + 1, size(C2, 2), size(C1, 2));
for i = 0 : n
    for j = 0 : n - i
        for a = 0 : n - i - j
            for b = 0 : n - i - j - a
                NSEC(1 + i, 1 + j, :, :) = NSEC(1 + i, 1 + j, :, :) + reshape(C1(1 + b + j, :)' * nchoosek(n, i) * nchoosek(n - i, j) * factorial(n - i - j) / (factorial(n - i - j - a - b) * factorial(a) * factorial(b)) * (- 1) ^ (a + b) * C2(1 + a + i, :), 1, 1, size(C2, 2), size(C1, 2));
            end
        end
    end
end
end