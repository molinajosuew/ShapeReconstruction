function GetSeparableBernsteinImageOut = NonSepBernImg(K, n, a, b)
I = zeros(a, b);
k = 1;
for i = 0 : n
    for j = 0 : n - i
        I = I + K(k) .* NonSepBern(1, a, b, n, i, j);
        k = k + 1;
    end
end
GetSeparableBernsteinImageOut = im2double(I <= 0);
end