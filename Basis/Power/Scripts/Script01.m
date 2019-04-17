% This script reconstructs one noisy power polynomial.

n = 4;

x_a = - 15;
x_b = 15;
y_a = - 15;
y_b = 15;

m_x = 100;
m_y = 100;

x_n = m_x * (x_b - x_a) + 1;
y_n = m_y * (y_b - y_a) + 1;

% P = GetRandomPolynomial(n, x_a, x_b, y_a, y_b, x_n, y_n, 10000, 10);

I = GetImageOfPowerPolynomial(P, x_a, x_b, y_a, y_b, x_n, y_n);
figure;
imshow(I);

[C, D] = PowerReconstruction(I, n, x_a, x_b, y_a, y_b, - 1);

figure;
imshow(D);

R = GetImageOfPowerPolynomial(C, x_a, x_b, y_a, y_b, x_n, y_n);
figure;
imshow(R);

figure;
imshow(abs(I - R));