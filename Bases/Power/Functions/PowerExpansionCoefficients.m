function PC = PowerExpansionCoefficients(n, a, b)
    N = n + 1;
    P = zeros(1, N + 1);
    for i = 0 : N
        P(1 + i) = nchoosek(N, i);
    end
    P = P * 2 ^ (- N + 1);
    PC = PolCoeffs(N, 0, P, a, b);
    for k = 1 : n
        PC(1 + k, :) = PolCoeffs(N, k, P, a, b);
    end
end