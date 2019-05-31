function SepAlgorithmOut = SeparableReconstruction(I, n)
    % System of SeparableMoment Equations
    A = zeros(2 * (n + 1) ^ 2, (n + 1) ^ 2);
    row = 1;
    for r = 0 : n
        for s = 0 : n
            for i = 0 : n
                for j = 0 : n
                    helper = SeparableProduct(n, i, r) * SeparableProduct(n, j, s);
                    A(row + 0, 1 + (n + 1) * i + j) = helper * (n + r) * (SeparableMoment(I, n + r - 1, n + s, i + r - 1, j + s) - SeparableMoment(I, n + r - 1, n + s, i + r, j + s));
                    A(row + 1, 1 + (n + 1) * i + j) = helper * (n + s) * (SeparableMoment(I, n + r, n + s - 1, i + r, j + s - 1) - SeparableMoment(I, n + r, n + s - 1, i + r, j + s));
                end
            end
            row = row + 2;
        end
    end
    % Power Basis Constraints
    B = zeros(1, size(A, 2));
    B(1, 1) = 1;
    row = 2;
    for r = 0 : 4
        for s = 0 : 4
            if r + s > 4
                B(row, :) = zeros(1, 25);
                for i = 0 : r
                    for j = 0 : s
                        B(row, 1 + (n + 1) * i + j) = nchoosek(4, r) * nchoosek(4, s) * nchoosek(r, i) * nchoosek(s, j) * (- 1) ^ (r + s - i - j);
                    end
                end
                row = row + 1;
            end
        end
    end
    SepAlgorithmOut = lsqlin(A, zeros(size(A, 1), 1), [], [], B, eye(size(B, 1), 1), [], [], []);
end