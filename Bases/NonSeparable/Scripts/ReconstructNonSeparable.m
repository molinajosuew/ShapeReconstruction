clear;
clc;

n = 4;
L = 2 ^ 5;
m_x = 2 ^ 4;
m_y = 2 ^ 4;
psnr = - 1;

x_n = m_x * L + 1;
y_n = m_y * L + 1;

[RC, D] = NonSeparableReconstruction(I, n, L, psnr);
R = GetImageOfNonSeparable(RC, L, x_n, y_n);

imshow(abs(I - R));