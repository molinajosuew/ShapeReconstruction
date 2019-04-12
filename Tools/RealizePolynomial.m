% INPUT
%     coeffs = coefficients of poly
% a, b, c, d = bounds of domain
%       r, s = resolution of domain
% OUTPUT
% img = r by s matrix taking values in {0, 1} of poly p^{-1}((-inf, 0]) \cap [a, b] × [c, d] induced by coeffs

function img = RealizePolynomial(C, xA, xB, yA, yB, xN, yN)
    [X, Y] = meshgrid(linspace(xA, xB, xN), linspace(yA, yB, yN));
    img = zeros(size(X));
    cnt = 1;
    n = (sqrt(8 * size(C, 2) + 1) - 3) / 2; % get degree of poly
    for i = 0 : n
        for j = 0 : n - i
            img = img + C(cnt) * X .^ i .* Y .^ j;
            cnt = cnt + 1;
        end
    end
    img = im2double(img <= 0)';
end