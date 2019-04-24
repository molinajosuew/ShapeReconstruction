function F = ProtoNonSeparableBernsteinFactor(n, a, b, i, j, L)
    F = L ^ (a + b) * factorial(n) * factorial(i + a) * factorial(j + b) / (factorial(i) * factorial(j) * factorial(n + a + b));
end