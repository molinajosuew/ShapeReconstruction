function I = RandomStablyBounded4thDegreeBivariatePolynomial(ra, rb, ca, cb, rN, cN, R)
    A = RandSymmPosDefMatrix(3);
    C = rand(1, 15) * R;
    C(1) = 1;
    [X, Y] = meshgrid(linspace(ca, cb, cN), linspace(ra, rb, rN));
    I = arrayfun(@(x, y) [x ^ 2, x * y, y ^ 2] * A * [x ^ 2, x * y, y ^ 2]', X, Y);
    cnt = 1;
    for i = 0 : 4
        for j = 0 : 4 - i
            if i + j < 4
                I = I + C(cnt) * X .^ i .* Y .^ j;
            else
                C(cnt) = 0;
            end
            cnt = cnt + 1;
        end
    end
    I = im2double(I <= 0);
end