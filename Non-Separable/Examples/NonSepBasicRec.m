clear;
clc;

n = 4;
load('Image01.mat');

K = NonSepAlg(ConvolveAndDownsample(I, 2 * n - 1, 1, 1), n);
R = NonSepBernImg(K, n, size(I, 1), size(I, 2));

imshow(cat(2, cat(2, cat(2, cat(2, I, repmat(255, 256, 10)), R), repmat(255, 256, 10)), abs(I - R)));