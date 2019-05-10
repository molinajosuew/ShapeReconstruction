function M = GetRandomSymmetricPositiveDefiniteMatrix()
    V = zeros(3, 3);

    V(1, :) = (2 * rand() - 1) * [1, 1, 1];
    V(2, :) = (2 * rand() - 1) * [1, 1, 0];
    V(3, :) = (2 * rand() - 1) * [1, 0, 0];

    M = zeros(3, 3);

    for i = 1 : 3
        for j = 1 : 3
            M(i, j) = dot(V(i, :), V(j, :));
        end
    end
end