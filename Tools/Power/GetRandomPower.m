function [P, I] = GetRandomPower(n, x_a, x_b, y_a, y_b, is_other_mode)
    [X, Y] = meshgrid(linspace(x_a, x_b, 2 ^ 9 + 1), linspace(y_a, y_b, 2 ^ 9 + 1));
    I = zeros(size(X));
    while I(1, 1) == 1 || length(unique(I(1, :))) ~= 1 || length(unique(I(end, :))) ~= 1 || length(unique(I(:, 1))) ~= 1 || length(unique(I(:, end))) ~= 1 || length(unique(I)) ~= 2
        I = zeros(size(X));
        P = zeros(1, nchoosek(n + 2, 2));
        k = 1;
        if ~ is_other_mode
            for i = 0 : n
                for j = 0 : n - i
                    P(k) = 2 * rand() - 1;
                    I = I + P(k) * X .^ i .* Y .^ j;
                    k = k + 1;
                end
            end
        else
            for l = 0 : n
                for i = 0 : l
                    j = l - i;
                    P(k) = 2 * rand() - 1;
                    I = I + P(k) * X .^ i .* Y .^ j;
                    k = k + 1;
                end
            end
        end
        I = im2double(I <= 0);
    end
end