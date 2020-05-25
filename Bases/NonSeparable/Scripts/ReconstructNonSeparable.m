clear;
clc;

stats = zeros(100, 50);

for image_n = 1 : 100
    disp(image_n);
    [~, I] = GetRandomNonSeparable(4, 3);
    
    for L = 1 : 50
        R = GetImageOfNonSeparable(NonSeparableReconstruction(I, 4, 1, - 1, L), 1, 257, 257, L);
        stats(image_n, L) = SorensenDiceCoefficient(I, R);
    end
end

boxplot(stats, 'whisker', 1000);
xlabel("L" + newline + "(as defined in (7) of our paper)");
ylabel('Sørensen–Dice Coefficient');
title('Bernstein Reconstruction of 100 Algebraic Shapes in [0,1]^2');