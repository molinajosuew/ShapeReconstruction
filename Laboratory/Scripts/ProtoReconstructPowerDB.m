n = 4;

x_a = - 175;
x_b = 175;
y_a = - 175;
y_b = 175;

m_x = 1;
m_y = 1;

psnr = - 1;

x_n = m_x * (x_b - x_a) + 1;
y_n = m_y * (y_b - y_a) + 1;

% P = GetRandomPowerPolynomial(n, x_a, x_b, y_a, y_b, false);
% I = GetImageOfPowerPolynomial(P, x_a, x_b, y_a, y_b, x_n, y_n);
figure;
imshow(I);

[C, D] = ProtoPowerReconstructionDB(I, n, x_a, x_b, y_a, y_b, psnr);

figure;
imshow(D);

R = GetImageOfPowerPolynomial(C, x_a, x_b, y_a, y_b, x_n, y_n);
figure;
imshow(R);

figure;
imshow(~abs(I - R));

disp(SorensenDiceCoefficient(I, R));