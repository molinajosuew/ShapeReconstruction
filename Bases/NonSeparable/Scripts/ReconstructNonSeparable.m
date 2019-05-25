clear;
clc;

n = 2;
L = 2;
m_x = 250;
m_y = 250;
psnr = - 1;

x_n = m_x * L + 1;
y_n = m_y * L + 1;
I = GetImageOfNonSeparable(GetRandomNonSeparable(n, L), L, x_n, y_n);
[C, D] = NonSeparableReconstruction(I, n, L, psnr);
R = GetImageOfNonSeparable(C, L, x_n, y_n);

figure;
imshow(D);
figure;
imshow(abs(I - R));