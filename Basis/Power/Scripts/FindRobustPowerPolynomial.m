% This script finds noisy images that reconstruct well.

clear;
clc;

n = 4;

x_a = - 5;
x_b = 5;
y_a = - 5;
y_b = 5;

m_x = 30;
m_y = 30;

psnr = - 1;

x_n = m_x * (x_b - x_a) + 1;
y_n = m_y * (y_b - y_a) + 1;

c = 1;

while c <= 10
    while true
        I = GetImageOfPowerPolynomial(GetRandomPowerPolynomial(n, x_a, x_b, y_a, y_b, false), x_a, x_b, y_a, y_b, x_n, y_n);
        R = GetImageOfPowerPolynomial(PowerReconstruction(I, n, x_a, x_b, y_a, y_b, psnr), x_a, x_b, y_a, y_b, x_n, y_n);
        if abs(1 - SorensenDicecoefficient(I, R)) < 0.01
            figure;
            imshow(abs(I - R));
            break;
        end
    end
    
    c = c + 1;
end