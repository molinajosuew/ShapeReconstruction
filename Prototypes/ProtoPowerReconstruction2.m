function R = ProtoPowerReconstruction2(I, n, mX, mY)
    J = conv2(wjmBSpline(n + ceil(n / 2), mY), wjmBSpline(n + ceil(n / 2), mX), I, "same");
    J = J(1 : mY : end, 1 : mX : end);

    C = ProtoPowerCoefficients(n + ceil(n / 2), - size(J, 1), size(J, 1));
    
    M = zeros(2 * (ceil(n / 2) + 1) ^ 2, nchoosek(n + 2, 2));
    
    row = 1;

    for r = 0 : ceil(n / 2)
        for s = 0 : ceil(n / 2)
            col = 1;

            for i = 0 : n
                for j = 0 : n - i
                    M(row + 0, col) = (i + r) * ProtoMoment(C, i + r - 1, j + s + 0, J);
                    M(row + 1, col) = (j + s) * ProtoMoment(C, i + r + 0, j + s - 1, J);

                    col = col + 1;
                end
            end

            row = row + 2;
        end
    end

    R = RealizePolynomial(lsqlin(M, zeros(size(M, 1), 1), [], [], eye(1, size(M, 2)), 1, [], [])', - size(J, 1), size(J, 1), - size(J, 1), size(J, 1), size(I, 2), size(I, 2));
%     R = RealizePolynomial(quadprog(M' * M, [], [], [], eye(1, size(M, 2)), 1)', 0, size(J, 1) - 1, 0, size(J, 1) - 1, size(I, 2), size(I, 2));
end