% n = degree
% a & b = resolution
function NSEC = NonSepExpCoeffs(n, a, b)
    K2 = DualExpansionCoefficients(n, a);
    K1 = DualExpansionCoefficients(n, b);
    NSEC = zeros(n + 1, n + 1, size(K1, 2), size(K2, 2));

    for i = 0 : n
        for j = 0 : n - i
            for q = 0 : n - i - j
                for r = 0 : n - i - j - q
                    NSEC(1 + i, 1 + j, :, :) = NSEC(1 + i, 1 + j, :, :) + reshape(K2(1 + r + j, :)' * nchoosek(n, i) * nchoosek(n - i, j) * factorial(n - i - j) / (factorial(n - i - j - q - r) * factorial(q) * factorial(r)) * (- 1) ^ (q + r) * K1(1 + q + i, :), 1, 1, size(K1, 2), size(K2, 2));
                end
            end
        end
    end
end