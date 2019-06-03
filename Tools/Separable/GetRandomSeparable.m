function [P, I] = GetRandomSeparable(n, x_a, x_b, y_a, y_b)
    [X, Y] = meshgrid(linspace(x_a, x_b, 2 ^ 9 + 1), linspace(y_a, y_b, 2 ^ 9 + 1));
    I = zeros(size(X));
    while I(1, 1) == 1 || length(unique(I(1, :))) ~= 1 || length(unique(I(end, :))) ~= 1 || length(unique(I(:, 1))) ~= 1 || length(unique(I(:, end))) ~= 1 || length(unique(I)) ~= 2
        I = zeros(size(X));
        P = zeros(1, (n + 1) ^ 2);
        k = 1;
        for i = 0 : n
            for j = 0 : n
                P(k) = 2 * rand() - 1;
                I = I + P(k) .* nchoosek(n, i) .* ((X - x_a) / (x_b - x_a)) .^ i .* (1 - (X - x_a) / (x_b - x_a)) .^ (n - i) .* nchoosek(n, j) .* ((Y - y_a) / (y_b - y_a)) .^ j .* (1 - (Y - y_a) / (y_b - y_a)) .^ (n - j);
                k = k + 1;
            end
        end
        I = im2double(I <= 0);
    end
end