function PowerBaseImgOut = PowerBaseImg(K, n, r, s, l)
I = zeros(r, s);
k = 1;
for i = 0 : n
    for j = 0 : n - i
        I = I + K(k) .* PowerBase(i, j, r, s, l);
        k = k + 1;
    end
end
PowerBaseImgOut = im2double(I <= 0);
end