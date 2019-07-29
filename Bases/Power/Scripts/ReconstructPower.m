clear;
clc;

n = 4;
x_a = - 1;
x_b = 1;
y_a = - 1;
y_b = 1;
m_x = 2 ^ 8;
m_y = 2 ^ 8;
psnr = - 1;

x_n = m_x * (x_b - x_a) + 1;
y_n = m_y * (y_b - y_a) + 1;

I = GetImageOfPower(GetRandomPower(n, x_a, x_b, y_a, y_b, false), x_a, x_b, y_a, y_b, x_n, y_n);
R = GetImageOfPower(PowerReconstruction(I, n, x_a, x_b, y_a, y_b, psnr), x_a, x_b, y_a, y_b, x_n, y_n);

imshow(abs(I - R));