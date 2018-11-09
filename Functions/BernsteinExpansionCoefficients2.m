function BEC2O = BernsteinExpansionCoefficients2(n, r)
DEC2O = DualExpansionCoefficients2(n, r);
BEC2O = zeros(n + 1, size(DEC2O, 2));
for v = 0 : n
    for k = 1 : size(DEC2O, 2)
        BEC2O(v + 1, k) = 0;
        for l = v : n
            BEC2O(v + 1, k) = BEC2O(v + 1, k) + nchoosek(n, l) * nchoosek(l, v) * (- 1) ^ (l - v) * DEC2O(l + 1, k);
        end
    end
end
end