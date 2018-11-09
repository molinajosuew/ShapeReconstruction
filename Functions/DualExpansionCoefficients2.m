function DEC2O = DualExpansionCoefficients2(n, r)
n = n + 1;
r = r - 1;
P = zeros(1, n + 1);
for i = 0 : n
    P(1 + i) = nchoosek(n, i);
end
P = P * 2 ^ (- n + 1);
DEC2O = zeros(n + 1, r + n - (mod(n, 2) == 0));
for k = 0 : n - 1
    DEC2O(1 + k, :) = PolCoeffs(n, k, P, 0, r) * r ^ (- k);
end
end