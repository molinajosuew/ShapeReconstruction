function N = ProtoNonSeparableBernsteinExpansionCoefficients(n, L)
    P = PowerExpansionCoefficients(n, - L, L);
    N = zeros(n, n, n, size(P, 2), size(P, 2));
    
    for k = 0 : n
        for i = 0 : k
            for j = 0 : k - i
                S = zeros(size(P, 2), size(P, 2));

                for a = 0 : k - i - j
                    for b = 0 : k - i - j - a
                        c = k - i - j - a - b;

                        S = S + nchoosek(k, i) .* nchoosek(k - i, j) .* factorial(k - i - j) .* L .^ (a - k) .* (- 1) .^ (b + c) .* (P(1 + c + j, :)' * P(1 + b + i, :)) .* (factorial(a) .* factorial(b) .* factorial(c)) .^ (- 1);
                    end
                end

                N(1 + k, 1 + i, 1 + j, :, :) = S;
            end
        end
    end
end