clear;
clc;

n = 4;
L = 20;
m_x = 25;
m_y = 25;
psnr = 50;
N = 100;

x_n = m_x * L + 1;
y_n = m_y * L + 1;
D = zeros(N, 4);

load('images01.mat');

for i = 1 : N
%     I = GetImageOfPowerPolynomial(GetRandomPowerPolynomial(n, 0, L, 0, L, false), 0, L, 0, L, x_n, y_n);
%     I = GetImageOfNonSeparableBernsteinPolynomial(GetRandomNonSeparableBernsteinPolynomial(n, L), L, x_n, y_n);
    I = reshape(II(i, :, :), [size(II, 2), size(II, 3)]);
    
    D(i, 1) = SorensenDiceCoefficient(I, GetImageOfPowerPolynomial(PowerReconstruction(I, n, 0, L, 0, L, psnr), 0, L, 0, L, x_n, y_n));
    
    D(i, 2) = SorensenDiceCoefficient(I, GetImageOfNonSeparableBernsteinPolynomial(ProtoNonSeparableBernsteinReconstruction(I, n, L, psnr), L, x_n, y_n));
    
    [QP_E, QP] = VetterliReconstruct(I, 0, L, 0, L, x_n, y_n, psnr);
    D(i, 3) = SorensenDiceCoefficient(I, QP);
    D(i, 4) = SorensenDiceCoefficient(I, QP_E);
    
    disp(100 * i / N + "%");
end

% figure;
% histogram(S(:, k));

figure;
boxplot(D);