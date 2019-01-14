function PowerBaseImgOut = PowerBaseImg(K, n, a, b)
I = zeros(a, b);
k = 1;
for i = 0 : n
    for j = 0 : n - i
        I = I + K(k) .* PowerBase(i, j, a, b);
        k = k + 1;
    end
end
PowerBaseImgOut = im2double(I <= 0);
end