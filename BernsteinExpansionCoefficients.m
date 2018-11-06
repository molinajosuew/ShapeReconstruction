function BernsteinExpansionCoefficientsOut = BernsteinExpansionCoefficients(n, r)
DualExpansionCoefficientsOut = DualExpansionCoefficients(n, r);
BernsteinExpansionCoefficientsOut = zeros(n + 1, size(DualExpansionCoefficientsOut, 2));
for v = 0 : n
    for k = 1 : size(DualExpansionCoefficientsOut, 2)
        BernsteinExpansionCoefficientsOut(v + 1, k) = 0;
        for l = v : n
            BernsteinExpansionCoefficientsOut(v + 1, k) = BernsteinExpansionCoefficientsOut(v + 1, k) + nchoosek(n, l) * nchoosek(l, v) * (- 1) ^ (l - v) * DualExpansionCoefficientsOut(l + 1, k);
        end
    end
end
end