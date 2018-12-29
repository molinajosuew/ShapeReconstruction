clear;
clc;

n = 4;
a = 128;
b = 128;

NSEC = NonSepExpCoeffs(n, a, b);

for i = 0 : n
    for j = 0 : n - i
        imshow(conv2(wjmBSpline(n, 1), wjmBSpline(n, 1), reshape(NSEC(1 + i, 1 + j, :, :), size(NSEC, 4), size(NSEC, 3)), 'valid'));
        pause(.5);
    end
end