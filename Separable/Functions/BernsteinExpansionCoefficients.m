function BECO = BernsteinExpansionCoefficients(n, r)
DEC2O = DualExpansionCoefficients(n, r, 1);
BECO = zeros(n + 1, size(DEC2O, 2));
for v = 0 : n
    for k = 1 : size(DEC2O, 2)
        BECO(v + 1, k) = 0;
        for l = v : n
            BECO(v + 1, k) = BECO(v + 1, k) + nchoosek(n, l) * nchoosek(l, v) * (- 1) ^ (l - v) * DEC2O(l + 1, k);
        end
    end
end
end