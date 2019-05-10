clear;
clc;

disp("getting polynomial");
n = 4;
k = 9;
[C1, a, b, c, d] = GetRandomBoundedPolynomial(n, 0, 1, 0, 1, 5 * k);
I = RealizePolynomial(C1, a, b, c, d, 2 ^ k + 1, 2 ^ k + 1);
figure;
imshow(I);

disp("down-sampling");
p = 4;
D = ConvolveAndDownsample(I, n + ceil(n / 2), 2 ^ p, 2 ^ p);
figure;
imshow(D / 2 ^ (2 * p));

disp("reconstructing");
C2 = PowerBaseAlg(D, n, 1);
R = RealizePolynomial(C2, 0, 1, 0, 1, 2 ^ k + 1, 2 ^ k + 1);
figure;
imshow(R);

figure;
imshow(abs(I - R));