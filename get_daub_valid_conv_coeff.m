function [C, F] = get_daub_valid_conv_coeff(n, a, b)
    n = n + 1;
    if mod(n, 2) == 0
        a = a - 1;
    end
    b = b + 1;
    P = 2 * dbwavf(['db', num2str(n)]);
    for k = 0 : n - 1
        C(k + 1, :) = PolCoeffs(n, k, P, a, b);
    end
    F = ScalingIntegers(P);
end