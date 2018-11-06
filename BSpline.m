% n - B-Spline Degree
% m - B-Spline Dilation Factor

function BSplineOut = BSpline(n, m)
X = floor(- (n * m + m) / 2) : ceil((n * m + m) / 2);
B1(1, :) = (- m / 2 < X & X < m / 2) * 1 + (- m / 2 == X | X == m / 2) * .5;
B2(1, :) = (- m < X & X < 0) * 1 + (- m == X | X == 0) * .5;
for i = 1 : n
    B1(i + 1, :) = ((X / m + (i + 1) / 2) .* B2(i, :) + ((i + 1) / 2 - X / m) .* [zeros(1, m), B2(i, 1 : end - m)]) / i;
    B2(i + 1, :) = ((X / m + (i + 2) / 2) .* [B1(i, m + 1 : end), zeros(1, m)] + (i / 2 - X / m) .* B1(i, :)) / i;
end
BSplineOut = B1(end, :);
BSplineOut = BSplineOut(2 : end - 1);
end