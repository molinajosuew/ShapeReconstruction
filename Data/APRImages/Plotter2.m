clear;
clc;

n = 4;
L = 30;
C = NonSeparableExpansionCoefficients(n, L);
C_c = (size(C, 4) + 1) / 2;

for i = 0 : n
    for j = 0 : n - i
        figure;
        I = reshape(C(1 + n, 1 + i, 1 + j, :, :), size(C, 4), size(C, 5));
        imagesc(I(C_c : C_c + 15, C_c : C_c + 15));
        
        xticks(1 : 16);
        xticklabels(0 : 15);
        yticks(1 : 16);
        yticklabels(0 : 15);
        
        set(gca,'YDir','normal');
        colorbar;
        title('Non-Separable Bernstein Expansion Coefficients \eta_{k, l}^{4, ' + string(i) + ', ' + string(j) + ', 4, 4, 30}');
        xlabel('k');
        ylabel('l');
    end
end