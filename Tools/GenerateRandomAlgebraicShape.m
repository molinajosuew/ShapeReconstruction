function I = GenerateRandomAlgebraicShape(n, a, b, c, d, r, s, S)
    [X, Y] = meshgrid(linspace(a, b, r), linspace(c, d, s));

    I = zeros(r, s);

    while I(1, 1) == 1 || length(unique(I(1, :))) ~= 1 || length(unique(I(end, :))) ~= 1 || length(unique(I(:, 1))) ~= 1 || length(unique(I(:, end))) ~= 1 || length(unique(I)) ~= 2
        I = zeros(r, s);
        coefficients = [];

        for i = 0 : n
            for j = 0 : n - i
                if i == 0 && j == 0
                    R = 1;
                else
                    R = randi([- S, S]);
                end
                coefficients = [coefficients, R];
                I = I + R * X .^ i .* Y .^ j;
            end
        end

        I = im2double(I <= 0);
    end
end