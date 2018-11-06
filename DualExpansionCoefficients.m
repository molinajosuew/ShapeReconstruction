% This function returns a matrix of n + 1 rows, each row k (zero-indexed)
% containing the dual coefficients required to reproduce the monomial x ^ k
% over the interval [0, 1] uniformly partitioned with r points.

function DualExpansionCoefficientsOut = DualExpansionCoefficients(n, r)
DualExpansionCoefficientsOut = StrangFix_coefficients(n, 1, [- floor(n / 2), r - 1 + floor(n / 2)], 0);
for k = 0 : n
    DualExpansionCoefficientsOut(k + 1, :) = DualExpansionCoefficientsOut(k + 1, :) * (r - 1) ^ (- k);
end
end