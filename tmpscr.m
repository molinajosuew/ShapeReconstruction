% Generate a ton of images that reconstruct well for us.
% Reconstruct with and record the errors of our method and their method.
% Boxplot that shit.

clear;
clc;

load('/home/wjm/Documents/GitHub/ShapeReconstruction/Data/Statistics/images01.mat');

n = 4;
L = 5;
x_a = - L;
x_b = L;
y_a = - L;
y_b = L;
m_x = 50;
m_y = 50;
psnr = 42.2;
x_n = m_x * (x_b - x_a) + 1;
y_n = m_y * (y_b - y_a) + 1;

for j = 1 : 100
    try
        I = reshape(II(1, :, :), size(II, 2), size(II, 3));
        % GET BEST!!!
%         R1 = GetImageOfPower(PowerReconstruction(I, n, x_a, x_b, y_a, y_b, psnr), x_a, x_b, y_a, y_b, x_n, y_n);
%         for i = 1 : 100
%             R1_t = GetImageOfPower(PowerReconstruction(I, n, x_a, x_b, y_a, y_b, psnr), x_a, x_b, y_a, y_b, x_n, y_n);
%             
%             if SorensenDiceCoefficient(I, R1_t) > SorensenDiceCoefficient(I, R1)
%                 R1 = R1_t;
%             end
%         end
        [~, R2, ~] = VetterliReconstruction(I, x_a, x_b, y_a, y_b, x_n, y_n, psnr);
%         D(j, 1) = SorensenDiceCoefficient(I, R1);
        DD(j, 2) = SorensenDiceCoefficient(I, R2);
    catch
        disp("There was an error with Vetterli.");
    end
    disp(j);
end

load('/home/wjm/Documents/GitHub/ShapeReconstruction/tmpdat.mat');
D = [D(:, 2), DD(:, 2)];
boxplot(D);