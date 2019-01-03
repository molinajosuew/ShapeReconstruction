function NonSepAlgorithmOut = NonSepAlg(I, n)
A = zeros(2 * (ceil(n / 2) + 1) ^ 2, nchoosek(n + 2, 2));
row = 1;
for r = 0 : ceil(n / 2)
    for s = 0 : ceil(n / 2)
        col = 1;
        for i = 0 : n
            for j = 0 : n - i
                helper_1 = n * NonSepBernProd(n, i, j, r, s);
                helper_2 = NonSepMoment(I, n + r + s - 1, i + r, j + s);
                A(row + 0, col) = helper_1 * (NonSepMoment(I, n + r + s - 1, i + r - 1, j + s) - helper_2);
                A(row + 1, col) = helper_1 * (NonSepMoment(I, n + r + s - 1, i + r, j + s - 1) - helper_2);
                col = col + 1;
            end
        end
        row = row + 2;
    end
end
B = zeros(1, size(A, 2));
B(1, 1) = 1;
NonSepAlgorithmOut = lsqlin(A, zeros(size(A, 1), 1), [], [], B, eye(size(B, 1), 1), [], []);
end