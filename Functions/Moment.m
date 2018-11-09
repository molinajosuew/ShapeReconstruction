function MomentOut = Moment(I, n, m, k, l)
if k >= 0 && l >= 0 && n >= k && m >= l
    K1 = BernsteinExpansionCoefficients2(n, size(I, 1));
    K2 = BernsteinExpansionCoefficients2(m, size(I, 2));
    K1 = K1(k + 1, 1 + floor(n / 2) : end - floor(n / 2));
    K2 = K2(l + 1, 1 + floor(m / 2) : end - floor(m / 2));
    MomentOut = sum(sum(K1' * K2 .* I));
else
    MomentOut = 0;
end
end