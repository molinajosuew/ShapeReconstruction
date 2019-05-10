function PC = ProtoPowerExpansionCoefficientsDB(n, a, b)
    N = n + 1;
    P = dbwavf("db" + n);
    PC = PolCoeffs(N, 0, P, a, b);
    for k = 1 : n
        PC(1 + k, :) = PolCoeffs(N, k, P, a, b);
    end
end