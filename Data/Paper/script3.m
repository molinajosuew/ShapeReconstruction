% Get Promising

clear;
clc;

n = 4;
L = 10;
x_a = - L;
x_b = L;
y_a = - L;
y_b = L;
m_x = 2 ^ 5;
m_y = 2 ^ 5;
psnr = 50;
x_n = m_x * (x_b - x_a) + 1;
y_n = m_y * (y_b - y_a) + 1;

while true
    IC = GetRandomPower(n, - 1, 1, - 1, 1, false);
    I = GetImageOfPower(IC, - 1, 1, - 1, 1, x_n, y_n);
    
    [RC_best, D] = PowerReconstruction(I, n, x_a, x_b, y_a, y_b, psnr);
    R_B_best = GetImageOfPower(RC_best, x_a, x_b, y_a, y_b, x_n, y_n);
    
    for i = 1 : 100
        RC = PowerReconstruction(I, n, x_a, x_b, y_a, y_b, psnr);
        R_B = GetImageOfPower(RC, x_a, x_b, y_a, y_b, x_n, y_n);
        
        if SorensenDiceCoefficient(I, R_B) > SorensenDiceCoefficient(I, R_B_best)
            RC_best = RC;
            R_B_best = R_B;
        elseif SorensenDiceCoefficient(I, ~ R_B) > SorensenDiceCoefficient(I, R_B_best)
            RC_best = RC;
            R_B_best = ~ R_B;
        end
    end
    
    [~, Vet_rec, ~] = VetterliReconstruction(I, x_a, x_b, y_a, y_b, x_n, y_n, psnr);
    
    if abs(SorensenDiceCoefficient(I, R_B_best) - SorensenDiceCoefficient(I, Vet_rec)) < 0.005
        break;
    end
    
    disp(SorensenDiceCoefficient(I, R_B_best));
    disp(SorensenDiceCoefficient(I, Vet_rec));
    disp("---");
end

figure;
imshow(abs(I - R_B_best));
disp(SorensenDiceCoefficient(I, R_B_best));
figure;
imshow(abs(I - Vet_rec));
disp(SorensenDiceCoefficient(I, Vet_rec));