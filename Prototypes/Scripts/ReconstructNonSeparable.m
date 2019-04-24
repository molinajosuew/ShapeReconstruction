% This script loads and reconstructs non-separable Bernstein polynomials.

clear;
clc;

n = 4;
L = 100;

m_x = 5;
m_y = 5;

psnr = - 1;

x_n = m_x * L + 1;
y_n = m_y * L + 1;

while true
%     I = GetImageOfPowerPolynomial(GetRandomPowerPolynomial(n, - 5, 5, - 5, 5, false), - 5, 5, - 5, 5, x_n, y_n);
    I = GetImageOfNonSeparableBernsteinPolynomial(GetRandomNonSeparableBernsteinPolynomial(n, L), L, x_n, y_n);
    R = GetImageOfNonSeparableBernsteinPolynomial(ProtoNonSeparableBernsteinReconstruction(I, n, L, psnr), L, x_n, y_n);
    
    if abs(1 - SorensenDiceCoefficient(I, R)) < 0.01
        imshow(abs(I - R));
        break;
    end
end