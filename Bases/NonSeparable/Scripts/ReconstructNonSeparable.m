clear;
clc;

n = 4;
L = 10;
m_x = 50;
m_y = 50;
psnr = - 1;

x_n = m_x * L + 1;
y_n = m_y * L + 1;
I = GetImageOfNonSeparable(GetRandomNonSeparable(n, L), L, x_n, y_n);
R = GetImageOfNonSeparable(NonSeparableReconstruction(I, n, L, psnr), L, x_n, y_n);

imshow(abs(I - R));