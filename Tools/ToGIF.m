% This script corroborates that the expansion coefficients work.

clear;
clc;

n = 4;
L = 2;
K = linspace(0, L, 512);

[XX, YY] = meshgrid(K, K);

k = 1;

for i = 0 : n
    for j = 0 : n - i
        f = figure;        
                I = nchoosek(n, i) .* nchoosek(n - i, j) .* XX .^ i .* YY .^ j .* (L - XX - YY) .^ (n - i - j) ./ L .^ n;
                imagesc(K, K, I .* (XX < flip(YY)));
                title(sprintf("B_{%i,%i}^{%i} Restricted to Triangle", i, j, n));        
%         I = nchoosek(n, i) .* nchoosek(n - i, j) .* XX .^ i .* YY .^ j .* (2 * L - XX - YY) .^ (n - i - j) ./ (2 * L) .^ n;
%         imagesc(K, K, I);
%         title(sprintf("B_{%i,%i}^{%i} Restricted to Square", i, j, n));
        
        ax = gca;
        ax.YDir = "normal";
        colorbar;
        caxis([0, 1]);
        xticks(0 : .2 : 2);
        yticks(0 : .2 : 2);
        
        [X, map] = rgb2ind(frame2im(getframe(f)), 256);        
        if k == 1
            imwrite(X, map, 'NonSepTrg.gif', 'gif', 'Loopcount', inf);
        else
            imwrite(X, map, 'NonSepTrg.gif', 'gif', 'WriteMode', 'append');
        end        
        k = k + 1;
    end
end