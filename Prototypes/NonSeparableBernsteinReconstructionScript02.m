n = 4;
L = 1;
m_x = 5;
m_y = 5;

x_n = m_x * L + 1;
y_n = m_y * L + 1;

I = GetImageOfPolynomial(P, 0, L, 0, L, x_n, y_n);
figure;
imshow(I);

[C, D] = ProtoNonSeparableBernsteinReconstruction(I, n, L);
figure;
imshow(D);

R = GetImageOfNonSeparableBernsteinPolynomial(C, L, x_n, y_n);
figure;
imshow(R);

figure;
imshow(abs(I - R));