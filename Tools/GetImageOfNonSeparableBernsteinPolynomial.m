function I = GetImageOfNonSeparableBernsteinPolynomial(C, L, x_n, y_n)
    [X, Y] = meshgrid(linspace(0, L, x_n), linspace(0, L, y_n));
    I = zeros(size(X));
    n = (sqrt(8 * size(C, 2) + 1) - 3) / 2; % deduce degree
    k = 1;
    for i = 0 : n
        for j = 0 : n - i
            I = I + C(k) * (2 * L) ^ (- n) * nchoosek(n, i) * nchoosek(n - i, j) * X .^ i .* Y .^ j .* (2 * L - X - Y) .^ (n - i - j);
            k = k + 1;
        end
    end
    I = im2double(I <= 0);
end