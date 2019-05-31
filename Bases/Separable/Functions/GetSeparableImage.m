function GetSeparableBernsteinImageOut = GetSeparableImage(C, n, m, r, s)
    X = linspace(0, 1, r);
    Y = linspace(0, 1, s);
    [YY, XX] = meshgrid(X, Y);
    B = zeros(r, s);
    k = 1;
    for i = 0 : n
        for j = 0 : m
            B = B + C(k) .* nchoosek(n, i) .* XX .^ i .* (1 - XX) .^ (n - i) .* nchoosek(m, j) .* YY .^ j .* (1 - YY) .^ (m - j);
            k = k + 1;
        end
    end
    GetSeparableBernsteinImageOut = im2double(B <= 0);
end