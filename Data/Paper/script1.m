clear;
clc;

load('/home/wjmolina/Documents/GitHub/ShapeReconstruction/Data/Paper/05/dat.mat', 'I');

psnr = 100;

for L = [1, 10, 20]
    % Power Reconstruction
    
    [C_1_best, D] = PowerReconstruction(I, 4, - L, L, - L, L, psnr);
    R_1_best = ~ GetImageOfPower(C_1_best, - L, L, - L, L, size(I, 1), size(I, 2));
    
    for i = 1 : 1000
        C_1 = PowerReconstruction(I, 4, - L, L, - L, L, psnr);
        R_1 = ~ GetImageOfPower(C_1, - L, L, - L, L, size(I, 1), size(I, 2));
        
        if SorensenDiceCoefficient(I, R_1) > SorensenDiceCoefficient(I, R_1_best)
            C_1_best = C_1;
            R_1_best = R_1;
        end
    end
    
    figure;
    imshow(abs(I - R_1_best));
    title("1, " + L);
    title("1, " + L);
    
    % Vetterli Reconstruction
    
    R_2 = VetterliReconstruction(I, - L, L, - L, L, size(I, 1), size(I, 2), psnr);
    figure;
    imshow(abs(I - R_2));
    title("2, " + L);
    title("2, " + L);
end

figure;
title("placeholder");