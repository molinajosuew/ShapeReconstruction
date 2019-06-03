clear;
clc;

n = 4;
L = 2 ^ 4;
m_x = 2 ^ 5;
m_y = 2 ^ 5;
psnr = - 1;

x_n = m_x * L + 1;
y_n = m_y * L + 1;
IC = GetRandomNonSeparable(n, L);
I = GetImageOfNonSeparable(IC, L, x_n, y_n);
while true
    [RC, D] = NonSeparableReconstruction(I, n, L, psnr);
    R = GetImageOfNonSeparable(RC, L, x_n, y_n);
    
    % figure;
    % imagesc(I);
    % title('Image');
    % colorbar;
    % axis off;
    %
    % figure;
    % imagesc(D);
    % title('Down-Sample');
    % colorbar;
    % axis off;
    
    %     figure;
    %     imagesc(abs(I - R));
    %     title('Non-Separable Reconstruction 50 dB');
    %     colorbar;
    %     axis off;
    if SorensenDiceCoefficient(I, R) >= .99
        break;
    end
end

% figure;
% imagesc(D);
% title('Down-Sample');
% colorbar;
% axis off;

figure;
imagesc(abs(I - R));
title('Non-Separable Reconstruction 30 dB');
colorbar;
axis off;