clear;
clc;

n = 100;
res = zeros(n, 8);

for j = 1 : n
    [~, I] = GetRandomNonSeparable(4, 3);
    for i = 7 : - 1 : 0
        l = 2 ^ i;
        R = GetImageOfNonSeparable(NonSeparableReconstruction(I, 4, l, - 1), 1, 257, 257);
        if SorensenDiceCoefficient(I, ~ R) > SorensenDiceCoefficient(I, R)
            R = ~ R;
        end
        res(j, i + 1) = SorensenDiceCoefficient(I, R);
    end
end

boxplot(res);
xlabel('Triangle (0, 0), (0, 2^{k-1}), (2^{k-1}, 0)');
ylabel('Sørensen–Dice Coefficient');