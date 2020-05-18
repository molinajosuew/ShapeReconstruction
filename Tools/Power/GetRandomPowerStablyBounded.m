function I = GetRandomPowerStablyBounded(x_a, x_b, y_a, y_b, x_n, y_n)
    A = GetRandomSymmetricPositiveDefiniteMatrix();
    [X, Y] = meshgrid(linspace(x_a, x_b, x_n), linspace(y_a, y_b, y_n));
    while true
        C = 2 * rand(1, 15) - 1;
        C(1) = 1;
        I = arrayfun(@(x, y) [x ^ 2, x * y, y ^ 2] * A * [x ^ 2, x * y, y ^ 2]', X, Y);
        cnt = 1;
        for i = 0 : 4
            for j = 0 : 4 - i
                if i + j < 4
                    I = I + C(cnt) * X .^ i .* Y .^ j;
                else
                    C(cnt) = 0;
                end
                cnt = cnt + 1;
            end
        end
        I = im2double(I <= 0);
        if ~ (I(1, 1) == 1 || length(unique(I(1, :))) ~= 1 || length(unique(I(end, :))) ~= 1 || length(unique(I(:, 1))) ~= 1 || length(unique(I(:, end))) ~= 1 || length(unique(I)) ~= 2)
            break;
        end
    end
end