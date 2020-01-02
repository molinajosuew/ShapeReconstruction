% clear;
% clc;

n = 4;
<<<<<<< HEAD
x_a = - 16;
x_b = 16;
y_a = - 16;
y_b = 16;
m_x = 2 ^ 4;
m_y = 2 ^ 4;
psnr = - 1;
=======
x_a = - 5;
x_b = 5;
y_a = - 5;
y_b = 5;
m_x = 2 ^ 6;
m_y = 2 ^ 6;
psnr = 80;
>>>>>>> 5ba0275c5fd06411a460cb8a8cf485b9d8d60e3f

x_n = m_x * (x_b - x_a) + 1;
y_n = m_y * (y_b - y_a) + 1;

<<<<<<< HEAD
% I = GetImageOfPower(IC, x_a, x_b, y_a, y_b, x_n, y_n);
[C, D] = PowerReconstruction(I, n, x_a, x_b, y_a, y_b, psnr);
R = GetImageOfPower(C, x_a, x_b, y_a, y_b, x_n, y_n);
=======
I = GetImageOfPower(IC, - 1, 1, - 1, 1, x_n, y_n);
>>>>>>> 5ba0275c5fd06411a460cb8a8cf485b9d8d60e3f

[V_best, ~, ~] = VetterliReconstruction(I, x_a, x_b, y_a, y_b, x_n, y_n, psnr);

imshow(abs(I - V_best));