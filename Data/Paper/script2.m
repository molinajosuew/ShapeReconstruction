clear;
clc;

n = 4;
x_a = - 10;
x_b = 10;
y_a = - 10;
y_b = 10;
m_x = 2 ^ 5;
m_y = 2 ^ 5;
psnr = 30;
x_n = m_x * (x_b - x_a) + 1;
y_n = m_y * (y_b - y_a) + 1;

load('/home/wjm/Documents/GitHub/ShapeReconstruction/Data/Paper/03/dat.mat', 'IC');
I = GetImageOfPower(IC, - 1, 1, - 1, 1, x_n, y_n);

[RC_best, D] = PowerReconstruction(I, n, x_a, x_b, y_a, y_b, psnr);
R_B_best = ~ GetImageOfPower(RC_best, x_a, x_b, y_a, y_b, x_n, y_n);

for i = 1 : 100    
    RC = PowerReconstruction(I, n, x_a, x_b, y_a, y_b, psnr);
    R_B = ~ GetImageOfPower(RC, x_a, x_b, y_a, y_b, x_n, y_n);
    
    if SorensenDiceCoefficient(I, R_B) < SorensenDiceCoefficient(I, R_B_best)
        RC_best = RC;
        R_B_best = R_B;
    end
end

figure;
imshow(abs(I - R_B_best));