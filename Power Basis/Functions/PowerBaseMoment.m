function PowerBaseMomentOut = PowerBaseMoment(i, j, I)
if 0 <= i && 0 <= j
    K1 = DualExpansionCoefficients(j, size(I, 2));
    K2 = DualExpansionCoefficients(i, size(I, 1));
    PowerBaseMomentOut = sum(sum(K1(1 + j, 1 + floor(j / 2) : end - floor(j / 2))' * K2(1 + i, 1 + floor(i / 2) : end - floor(i / 2)) .* I));
else
    PowerBaseMomentOut = 0;
end
end