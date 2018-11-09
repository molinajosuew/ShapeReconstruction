function BECO = BernsteinExpansionCoefficients(n, r)
DECO = DualExpansionCoefficients(n, r);
BECO = zeros(n + 1, size(DECO, 2));
for v = 0 : n
    for k = 1 : size(DECO, 2)
        BECO(v + 1, k) = 0;
        for l = v : n
            BECO(v + 1, k) = BECO(v + 1, k) + nchoosek(n, l) * nchoosek(l, v) * (- 1) ^ (l - v) * DECO(l + 1, k);
        end
    end
end
end