function [P, I] = GetRandomPolynomial(n, x_a, x_b, y_a, y_b, x_n, y_n, r, pad)
    [X, Y] = meshgrid(linspace(x_a, x_b, x_n), linspace(y_a, y_b, y_n));
    I = zeros(size(X));
    while I(1, 1) == 1 || length(unique(I(1 : pad, :))) ~= 1 || length(unique(I(end + 1 - pad : end, :))) ~= 1 || length(unique(I(:, 1 : pad))) ~= 1 || length(unique(I(:, end + 1 - pad : end))) ~= 1 || length(unique(I)) ~= 2
        I = zeros(size(X));
        P = zeros(1, nchoosek(n + 2, 2));
        k = 1;
        for i = 0 : n
            for j = 0 : n - i
                P(k) = (2 * rand() - 1) * r;
                I = I + P(k) * X .^ i .* Y .^ j;
                k = k + 1;
            end
        end
        I = im2double(I <= 0);
    end
end