function PC = ProtoPowerExpansionCoefficientsDB(n, a, b)
    N = n + 1;
    W = 2 * dbaux(N);
    PC = PolCoeffs(N, 0, W, a, b);
    for k = 1 : n
        PC(1 + k, :) = PolCoeffs(N, k, W, a, b);
    end
end