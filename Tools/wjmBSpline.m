function wjmBSplineOut = wjmBSpline(n, m)
    X = - ceil((n * m + m) / 2) : ceil((n * m + m) / 2);
    B1(1, :) = (- m / 2 < X & X < m / 2) * 1 + (- m / 2 == X | X == m / 2) * .5;
    B2(1, :) = (- m < X & X < 0) * 1 + (- m == X | X == 0) * .5;
    for i = 1 : n
        B1(i + 1, :) = ((X / m + (i + 1) / 2) .* B2(i, :) + ((i + 1) / 2 - X / m) .* [zeros(1, m), B2(i, 1 : end - m)]) / i;
        B2(i + 1, :) = ((X / m + (i + 2) / 2) .* [B1(i, m + 1 : end), zeros(1, m)] + (i / 2 - X / m) .* B1(i, :)) / i;
    end
    wjmBSplineOut = B1(end, :);
    wjmBSplineOut = wjmBSplineOut(2 : end - 1);
end