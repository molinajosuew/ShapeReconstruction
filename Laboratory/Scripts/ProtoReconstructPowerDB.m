clear;
clc;

n = 4;

x_a = 0;
x_b = 2 ^ 9;
y_a = 0;
y_b = 2 ^ 9;

m_x = 2 ^ 0;
m_y = 2 ^ 0;

psnr = - 1;

x_n = m_x * (x_b - x_a) + 1;
y_n = m_y * (y_b - y_a) + 1;

P = GetRandomNonSeparable(n, x_b);
I = GetImageOfNonSeparable(P, x_b, x_n, y_n);

figure;
for t = 591 : 591
    [C, D] = ProtoPowerReconstructionDB(I, n, psnr, t);
    R = GetImageOfPower(C, x_a, x_b, y_a, y_b, x_n, y_n);
    imshow(abs(I - R));
    disp(t);
    pause(1);
end

figure;
for t = 595 : 595
    [C, D] = ProtoPowerReconstructionDB(I, n, psnr, t);
    R = GetImageOfPower(C, x_a, x_b, y_a, y_b, x_n, y_n);
    imshow(abs(I - R));
    disp(t);
    pause(1);
end