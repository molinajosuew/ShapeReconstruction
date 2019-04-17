function I = GetImageOfPowerPolynomial(C, x_a, x_b, y_a, y_b, x_n, y_n)
    [X, Y] = meshgrid(linspace(x_a, x_b, x_n), linspace(y_a, y_b, y_n));
    I = zeros(size(X));
    n = (sqrt(8 * size(C, 2) + 1) - 3) / 2; % deduce degree of poly
    k = 1;
    for i = 0 : n
        for j = 0 : n - i
            I = I + C(k) * X .^ i .* Y .^ j;
            k = k + 1;
        end
    end
    I = im2double(I <= 0);
end