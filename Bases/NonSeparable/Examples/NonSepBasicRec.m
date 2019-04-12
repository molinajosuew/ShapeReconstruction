clear;
clc;

disp("getting polynomial");
n = 4;
k = 8;
[C1, a, b, c, d] = GetRandomBoundedPolynomial(n, 0, .5, 0, .5, 2 * k);
I = RealizePolynomial(C1, a, b, c, d, 2 ^ k + 1, 2 ^ k + 1);
figure;
imshow(I);

disp("down-sampling");
p = 0;
D = ConvolveAndDownsample(I, 2 * n - 1, 2 ^ p, 2 ^ p);
figure;
imshow(D / 2 ^ (2 * p));

disp("reconstructing");
C2 = NonSepAlg(D, n);
R = NonSepBernImg(C2, n, 2 ^ k + 1, 2 ^ k + 1);
figure;
imshow(R);

figure;
imshow(abs(I - R));