function C = ProtoPowerExpansionCoefficientsDB(n, a, b)
    P = 2 * dbwavf("db" + num2str(n));
    for k = 0 : n - 1
        C(1 + k, :) = PolCoeffs(n, k, P, a - 1, b + 1 - mod(n, 2));
    end
end