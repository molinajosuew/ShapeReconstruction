clear;
clc;

n = 1000;
res = zeros(n, 9);

for j = 1 : n
    [~, I] = GetRandomPower(4, - 2, 2, - 2, 2, false);
    for i = 0 : 8
        l = 2 ^ i;
        R = GetImageOfPower(PowerReconstruction(I, 4, - l, l, - l, l, - 1), - l, l, - l, l, 513, 513);
        if SorensenDiceCoefficient(I, ~ R) > SorensenDiceCoefficient(I, R)
            R = ~ R;
        end
        res(j, i + 1) = SorensenDiceCoefficient(I, R);
    end
end

boxplot(res);
xlabel('Triangle (0, 0), (0, 2^{k-1}), (2^{k-1}, 0)');
ylabel('Sørensen–Dice Coefficient');