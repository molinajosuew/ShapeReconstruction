function BernsteinProductOut = BernsteinProduct(n, i, r)
BernsteinProductOut = 1;
if r > 0
    for k = 1 : r
        BernsteinProductOut = BernsteinProductOut * (i + k) / (n + k);
    end
end
end