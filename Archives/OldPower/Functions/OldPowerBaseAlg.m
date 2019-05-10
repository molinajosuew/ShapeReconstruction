function PowerBaseAlgOut = PowerBaseAlg(J, n, l)
A = zeros(2 * (ceil(n / 2) + 1) ^ 2, nchoosek(n + 2, 2));
row = 1;
for r = 0 : ceil(n / 2)
    for s = 0 : ceil(n / 2)
        col = 1;
        for i = 0 : n
            for j = 0 : n - i
                A(row + 0, col) = (i + r) * PowerBaseMoment(n + ceil(n / 2), i + r - 1, j + s, J, l);
                A(row + 1, col) = (j + s) * PowerBaseMoment(n + ceil(n / 2), i + r, j + s - 1, J, l);
                col = col + 1;
            end
        end
        row = row + 2;
    end
end
PowerBaseAlgOut = transpose(lsqlin(A, zeros(size(A, 1), 1), [], [], eye(1, size(A, 2)), 1, [], []));
% PowerBaseAlgOut = transpose(quadprog(transpose(A) * A, [], [], [], eye(1, size(A, 2)), 1));
end