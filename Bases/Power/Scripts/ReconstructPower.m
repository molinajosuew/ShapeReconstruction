% clear;
% clc;

n = 4;
x_a = - 5;
x_b = 5;
y_a = - 5;
y_b = 5;
m_x = 2 ^ 6;
m_y = 2 ^ 6;
psnr = 80;

x_n = m_x * (x_b - x_a) + 1;
y_n = m_y * (y_b - y_a) + 1;

I = GetImageOfPower(IC, - 1, 1, - 1, 1, x_n, y_n);

[V_best, ~, ~] = VetterliReconstruction(I, x_a, x_b, y_a, y_b, x_n, y_n, psnr);

imshow(abs(I - V_best));