function NonSepMomentOut = NonSepMoment(I, n, i, j)
if 0 <= n && 0 <= i && 0 <= j && i + j <= n
    K = NonSepExpCoeffs(n, size(I, 1), size(I, 2));
    NonSepMomentOut = sum(sum(reshape(K(1 + i, 1 + j, 1 + floor(n / 2) : end - floor(n / 2), 1 + floor(n / 2) : end - floor(n / 2)), size(I, 1), size(I, 2)) .* I));
else
    NonSepMomentOut = 0;
end
end