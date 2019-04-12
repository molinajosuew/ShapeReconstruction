function [R, D2] = ProtoPowerReconstruction(I, n, ra, rb, ca, cb)
    rN = size(I, 1);
    cN = size(I, 2);

    Tr = (rN - 1) / (rb - ra);
    Tc = (cN - 1) / (cb - ca);

    if floor(Tc) ~= Tc || floor(Tr) ~= Tr
        disp("image domain and resolution inconsistency"); % This can be fixed by up-sampling and shifting splines.
        return;
    end

    D1 = conv2(wjmBSpline(n, Tc), wjmBSpline(n, Tr), I) / (Tr * Tc);
    D1Cr = rN + (size(D1, 1) - rN) / 2;
    D1Cc = (size(D1, 2) - cN + 2) / 2;

    D2 = D1(D1Cr + (size(D1, 1) - rN) / 2 : - Tr : 1, D1Cc - Tc * floor((D1Cc - 1) / Tc) : Tc : end);

    C = ProtoPowerCoefficients(n + ceil(n / 2), - max(abs([ra, rb, ca, cb])), max(abs([ra, rb, ca, cb])));
    CC = (size(C, 2) + 1) / 2;

    M = zeros(2 * (ceil(n / 2) + 1) ^ 2, nchoosek(n + 2, 2));
    row = 1;

    for r = 0 : ceil(n / 2)
        for s = 0 : ceil(n / 2)
            col = 1;

            for i = 0 : n
                for j = 0 : n - i
                    AR = CC + (ra - rb) / (rN - 1) * (size(D2, 1) - D1Cr) + ra : (ra - rb) / (rN - 1) * (1 - D1Cr) + ra;
                    AC = CC + ca - floor((D1Cc - 1) / Tc) : CC + size(D2, 2) + ca - floor((D1Cc - 1) / Tc) - 1;

                    if 0 <= i + r - 1 && 0 <= j + s + 0
                        M(row + 0, col) = (i + r) * sum(sum(flip(C(i + r + 0, AR))' * C(j + s + 1, AC) .* D2));
                    else
                        M(row + 0, col) = 0;
                    end

                    if 0 <= i + r + 0 && 0 <= j + s - 1
                        M(row + 1, col) = (j + s) * sum(sum(flip(C(i + r + 1, AR))' * C(j + s + 0, AC) .* D2));
                    else
                        M(row + 1, col) = 0;
                    end

                    col = col + 1;
                end
            end

            row = row + 2;
        end
    end

    R = quadprog(M' * M, [], [], [], eye(1, size(M, 2)), 1)';
%     R = lsqlin(M, zeros(size(M, 1), 1), [], [], eye(1, size(M, 2)), 1, [], [])';
end