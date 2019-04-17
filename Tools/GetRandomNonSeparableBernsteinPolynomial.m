function [P, I] = GetRandomNonSeparableBernsteinPolynomial(n, L)
    [X, Y] = meshgrid(linspace(0, L, 2 ^ 7 + 1), linspace(0, L, 2 ^ 7 + 1));
    I = zeros(size(X));
    while I(1, 1) == 1 || length(unique(I(1, :))) ~= 1 || length(unique(I(end, :))) ~= 1 || length(unique(I(:, 1))) ~= 1 || length(unique(I(:, end))) ~= 1 || length(unique(I)) ~= 2
        I = zeros(size(X));
        P = zeros(1, nchoosek(n + 2, 2));
        k = 1;
        for i = 0 : n
            for j = 0 : n - i
                P(k) = 2 * rand() - 1;
                I = I + P(k) * nchoosek(n, i) * nchoosek(n - i, j) * X .^ i .* Y .^ j .* (2 * L - X - Y) .^ (n - i - j) / (2 * L) ^ n;
                k = k + 1;
            end
        end
        I = im2double(I <= 0);
    end
end