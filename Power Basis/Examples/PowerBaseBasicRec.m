clear;
clc;

n = 4;
l = 1;
load('Image01.mat');

K = PowerBaseAlg(ConvolveAndDownsample(I, n + ceil(n / 2), 1, 1), n, l);
R = PowerBaseImg(K, n, size(I, 1), size(I, 2), l);

imshow(cat(2, cat(2, cat(2, cat(2, I, repmat(255, 256, 10)), R), repmat(255, 256, 10)), abs(I - R)));