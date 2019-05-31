function SEC = SeparableExpansionCoefficients(n, a, b)
    PEC = PowerExpansionCoefficients(n, a, b);
    SEC = zeros(1 + n, size(PEC, 2));
    for k = 0 : n
        for r = 0 : k
            for s = 0 : n - k
                SEC(1 + k, :) = SEC(1 + k, :) + nchoosek(n, k) * nchoosek(k, r) * nchoosek(n - k, s) * (b - a) ^ (- n) * a ^ (k - r) * b ^ (n - k - s) * (- 1) ^ (k + r + s) * PEC(1 + r + s, :);
            end
        end
    end
end