clear;
clc;

n = 4;
L = 20;
m_x = 25;
m_y = 25;
psnr = 50;

x_n = m_x * L + 1;
y_n = m_y * L + 1;
IC = GetRandomNonSeparable(n, L);
I = GetImageOfNonSeparable(IC, L, x_n, y_n);
[RC, D] = NonSeparableReconstruction(I, n, L, psnr);
R = abs(I - GetImageOfNonSeparable(RC, L, x_n, y_n));

% figure;
% imshow(D);
figure;
imshow(R);