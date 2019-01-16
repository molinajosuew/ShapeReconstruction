function I = RandomStablyBounded4thDegreeBivariatePolynomial(a, b, c, d, k, l)
A = RandSymmPosDefMatrix(3);
C = (2 * rand(1, 15) - 1) * 100;
C(1) = 1;
[XX, YY] = meshgrid(linspace(a, b, k), linspace(c, d, l));
I = arrayfun(@(x, y) [x ^ 2, x * y, y ^ 2] * A * [x ^ 2, x * y, y ^ 2]', XX, YY);
k = 1;
for i = 0 : 4
    for j = 0 : 4 - i
        if i + j < 4
            I = I + C(k) * XX .^ i .* YY .^ j;
        else
            C(k) = 0;
        end
        k = k + 1;
    end
end
I = im2double(I <= 0);
end