% clear;
% clc;

n = 4;
x_a = - 16;
x_b = 16;
y_a = - 16;
y_b = 16;
m_x = 2 ^ 4;
m_y = 2 ^ 4;
psnr = - 1;

x_n = m_x * (x_b - x_a) + 1;
y_n = m_y * (y_b - y_a) + 1;

% I = GetImageOfPower(IC, x_a, x_b, y_a, y_b, x_n, y_n);
[C, D] = PowerReconstruction(I, n, x_a, x_b, y_a, y_b, psnr);
R = GetImageOfPower(C, x_a, x_b, y_a, y_b, x_n, y_n);

imshow(abs(I - ~ R));