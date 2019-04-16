function N = ProtoNonSeparableBernsteinExpansionCoefficients(n, L)
    P = PowerExpansionCoefficients(n, 0, L);
    N = zeros(n, n, size(P, 2), size(P, 2));

    for i = 0 : n
        for j = 0 : n - i
            S = zeros(size(P, 2), size(P, 2));

            for a = 0 : n - i - j
                for b = 0 : n - i - j - a
                    c = n - i - j - a - b;

                    S = S + nchoosek(n, i) .* nchoosek(n - i, j) .* factorial(n - i - j) .* L .^ (a - n) .* (- 1) .^ (b + c) .* (P(1 + c + j, :)' * P(1 + b + i, :)) .* (factorial(a) .* factorial(b) .* factorial(c)) .^ (- 1);
                end
            end

            N(1 + i, 1 + j, :, :) = S;
        end
    end
end