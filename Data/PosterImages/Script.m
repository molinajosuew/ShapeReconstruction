clear;
clc;

n = 4;
L = 10;
x_n = 2 ^ 9 + 1;
y_n = 2 ^ 9 + 1;

C = GetRandomNonSeparable(n, L);
I = GetImageOfNonSeparable(C, L, x_n, y_n);

[R_E, R, D] = VetterliReconstruct(I, 0, L, 0, L, x_n, y_n, 50);

figure;
imagesc(D);
colorbar;
axis off;

figure;
imagesc(abs(I - R_E));
axis off;

% figure;
% imagesc(abs(I - R_E));
% axis off;