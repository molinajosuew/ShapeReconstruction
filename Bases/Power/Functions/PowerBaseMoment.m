function PowerBaseMomentOut = PowerBaseMoment(n, i, j, J, l)
if 0 <= i && 0 <= j
    K1 = DualExpansionCoefficients(n, size(J, 1), l);
    K2 = DualExpansionCoefficients(n, size(J, 2), l);
    PowerBaseMomentOut = sum(sum(K1(1 + j, 1 + floor(n / 2) : end - floor(n / 2))' * K2(1 + i, 1 + floor(n / 2) : end - floor(n / 2)) .* J));
else
    PowerBaseMomentOut = 0;
end
end