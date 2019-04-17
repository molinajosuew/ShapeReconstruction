% This script finds noisy images that reconstruct well.

clear;
clc;

n = 4;

x_a = - 250;
x_b = 250;
y_a = - 250;
y_b = 250;

m_x = 1;
m_y = 1;

x_n = m_x * (x_b - x_a) + 1;
y_n = m_y * (y_b - y_a) + 1;

c = 1;

while c <= 10
    while true
        I = GetImageOfPolynomial(GetRandomPolynomial(n, x_a, x_b, y_a, y_b, x_n, y_n, 1000, 10), x_a, x_b, y_a, y_b, x_n, y_n);
        E = abs(I - GetImageOfPolynomial(PowerReconstruction(I, n, x_a, x_b, y_a, y_b, 100), x_a, x_b, y_a, y_b, x_n, y_n));
        if sum(sum(E)) / size(E, 1) / size(E, 2) < 0.01
            figure;
            imshow(E);
            break;
        end
    end
    
    c = c + 1;
end