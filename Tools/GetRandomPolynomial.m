%          n = max degree of poly
% a, b, c, d = bounds of domain
%        pad = clearance of image with border

function C = GetRandomPolynomial(n, x_a, x_b, y_a, y_b, x_n, y_n, r, pad)
    [X, Y] = meshgrid(linspace(x_a, x_b, x_n), linspace(y_a, y_b, y_n));
    img = zeros(size(X));
    while img(1, 1) == 1 || length(unique(img(1 : pad, :))) ~= 1 || length(unique(img(end + 1 - pad : end, :))) ~= 1 || length(unique(img(:, 1 : pad))) ~= 1 || length(unique(img(:, end + 1 - pad : end))) ~= 1 || length(unique(img)) ~= 2
        img = zeros(size(X));
        C = zeros(1, nchoosek(n + 2, 2));
        k = 1;
        for i = 0 : n
            for j = 0 : n - i
                C(k) = (2 * rand() - 1) * r;
                img = img + C(k) * X .^ i .* Y .^ j;
                k = k + 1;
            end
        end
        img = im2double(img <= 0);
    end
end