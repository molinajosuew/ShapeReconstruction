% This script randomly generates and saves bounded non-separable Bernstein polynomials.

clear;
clc;

% Parameters

n = 4;
L = 100;
m_x = 5;
m_y = 5;

x_n = m_x * L + 1;
y_n = m_y * L + 1;

c = 1;
d = 9;

while c <= d
    while true
        % Image of Shape to Reconstruct

        [P, I] = GetRandomNonSeparableBernsteinPolynomial(n, L, x_n, y_n, 1000, 15);

        % Reconstruction and Image of Down-Sample

        [C, D] = ProtoNonSeparableBernsteinReconstruction(I, n, L);

        % Image of Reconstruction

        R = GetImageOfNonSeparableBernsteinPolynomial(C, L, x_n, y_n);

        % Error

        E = abs(I - R);
        
        if sum(sum(E)) / size(E, 1) / size(E, 2) < 0.01
            save("NSBP0" + c, 'P');
            break;
        end
    end
    
    c = c + 1;
end