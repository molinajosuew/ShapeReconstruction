function BernsteinFactorOut = SeparableFactor(n, i, r)
    BernsteinFactorOut = 1;
    if r > 0
        for k = 1 : r
            BernsteinFactorOut = BernsteinFactorOut * (i + k) / (n + k);
        end
    end
end