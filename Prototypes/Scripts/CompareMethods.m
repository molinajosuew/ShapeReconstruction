clear;
clc;

n = 4;
L = 50;
m_x = 50;
m_y = 50;
psnr = 50;
N = 100;

x_n = m_x * L + 1;
y_n = m_y * L + 1;
S = zeros(N, 4);

for i = 1 : N
%     I = GetImageOfPowerPolynomial(GetRandomPowerPolynomial(n, 0, L, 0, L, false), 0, L, 0, L, x_n, y_n);    
    I = GetImageOfNonSeparableBernsteinPolynomial(GetRandomNonSeparableBernsteinPolynomial(n, L), L, x_n, y_n);
    
    S(i, 1) = SorensenDiceCoefficient(I, GetImageOfPowerPolynomial(PowerReconstruction(I, n, 0, L, 0, L, psnr), 0, L, 0, L, x_n, y_n));
    
    S(i, 2) = SorensenDiceCoefficient(I, GetImageOfNonSeparableBernsteinPolynomial(ProtoNonSeparableBernsteinReconstruction(I, n, L, psnr), L, x_n, y_n));
    
    [QP_E, QP] = VetterliReconstruct(I, 0, L, 0, L, x_n, y_n, psnr);
    S(i, 3) = SorensenDiceCoefficient(I, QP);
    S(i, 4) = SorensenDiceCoefficient(I, QP_E);
end

figure;
histogram(S(:, 1));
figure;
histogram(S(:, 2));
figure;
histogram(S(:, 3));
figure;
histogram(S(:, 4));

figure;
boxplot(S);