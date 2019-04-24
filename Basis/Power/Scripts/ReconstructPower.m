% This script reconstructs one noisy power polynomial.

n = 4;

x_a = - 5;
x_b = 5;
y_a = - 5;
y_b = 5;

m_x = 25;
m_y = 30;

psnr = - 1;

x_n = m_x * (x_b - x_a) + 1;
y_n = m_y * (y_b - y_a) + 1;

% I = RandomStablyBounded4thDegreeBivariatePolynomial(x_a, x_b, y_a, y_b, x_n, y_n);

P = GetRandomPowerPolynomial(n, x_a, x_b, y_a, y_b, false);
I = GetImageOfPowerPolynomial(P, x_a, x_b, y_a, y_b, x_n, y_n);
% figure;
% imshow(I);

[C, D] = PowerReconstruction(I, n, x_a, x_b, y_a, y_b, psnr);

% figure;
% imshow(D);

R = GetImageOfPowerPolynomial(C, x_a, x_b, y_a, y_b, x_n, y_n);
% figure;
% imshow(R);

figure;
imshow(abs(I - R));