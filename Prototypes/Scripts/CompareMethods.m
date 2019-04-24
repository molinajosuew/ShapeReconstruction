clear;
clc;

n = 4;
L = 5;
m_x = 108;
m_y = 108;
psnr = 50;
N = 100;

x_n = m_x * L + 1;
y_n = m_y * L + 1;
S = zeros(N, 2);

for i = 1 : N
%     C = GetRandomPowerPolynomial(n, 0, L, 0, L, false);
%     I = GetImageOfPowerPolynomial(C, 0, L, 0, L, x_n, y_n);
    
    C = GetRandomNonSeparableBernsteinPolynomial(n, L);
    I = GetImageOfNonSeparableBernsteinPolynomial(C, L, x_n, y_n);
    
    S(i, 1) = SorensenDicecoefficient(I, GetImageOfPowerPolynomial(PowerReconstruction(I, n, 0, L, 0, L, psnr), 0, L, 0, L, x_n, y_n));
    S(i, 2) = SorensenDicecoefficient(I, GetImageOfNonSeparableBernsteinPolynomial(ProtoNonSeparableBernsteinReconstruction(I, n, L, psnr), L, x_n, y_n));
end

figure;
hist(S(:, 1));
figure;
hist(S(:, 2));