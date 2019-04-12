% INPUT
%     coeffs = coefficients of poly
% a, b, c, d = bounds of domain
%       r, s = resolution of domain
% OUTPUT
% img = r by s matrix taking values in {0, 1} of poly p^{-1}((-inf, 0]) \cap [a, b] × [c, d] induced by coeffs

function img = RealizePolynomial(C, x_a, x_b, y_a, y_b, x_n, y_n)
    [X, Y] = meshgrid(linspace(x_a, x_b, x_n), linspace(y_a, y_b, y_n));
    img = zeros(size(X));
    n = (sqrt(8 * size(C, 2) + 1) - 3) / 2; % deduce degree of poly
    k = 1;
    for i = 0 : n
        for j = 0 : n - i
            img = img + C(k) * X .^ i .* Y .^ j;
            k = k + 1;
        end
    end
    img = im2double(img <= 0);
end