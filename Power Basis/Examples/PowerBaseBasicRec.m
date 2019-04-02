clear;
clc;

l = 1;
n = 4;
load('Image01.mat');

% Images should have dimensions of the form 2 ^ k + 1, k > 0, to make the down-sampling work.
I = [I, zeros(size(I, 1), 1)];
I = [I; zeros(1, size(I, 2))];

K = PowerBaseAlg(ConvolveAndDownsample(I, n + ceil(n / 2), 16, 16), n, l);
J = PowerBaseImg(K, n, size(I, 1), size(I, 2), l);

figure;
imshow(I);
figure;
imshow(J);
figure;
imshow(abs(I - J));