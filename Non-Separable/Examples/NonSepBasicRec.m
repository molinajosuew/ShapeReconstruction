clear;
clc;

n = 4;
load('Image01.mat');

K = NonSepAlg(ConvolveAndDownsample(I, 2 * n - 1, 1, 1), n);
J = NonSepBernImg(K, n, size(I, 1), size(I, 2));

figure;
imshow(I);
figure;
imshow(J);
figure;
imshow(abs(I - J));