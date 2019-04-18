% This script loads and reconstructs non-separable Bernstein polynomials.

clear;
clc;

% Parameters

n = 4;
L = 2 ^ 5 + 1 - 5;
m_x = 108;
m_y = 108;

x_n = m_x * L + 1;
y_n = m_y * L + 1;

% Image of Shape to Reconstruct

P = GetRandomNonSeparableBernsteinPolynomial(n, L);

I = GetImageOfNonSeparableBernsteinPolynomial(P, L, x_n, y_n);
figure;
imshow(I);

% Reconstruction and Image of Down-Sample

[C, D] = ProtoNonSeparableBernsteinReconstruction(I, n, L, 100);
figure;
imshow(D);

% Error of Reconstruction

E = abs(I - GetImageOfNonSeparableBernsteinPolynomial(C, L, x_n, y_n));
figure;
imshow(E);