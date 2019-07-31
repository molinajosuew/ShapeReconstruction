function [C, D] = ProtoPowerReconstructionDB(I, n, psnr, t)
    S = ScalingIntegers(2 * dbwavf("db" + (n + ceil(n / 2) + 1)));
    D = conv2(S, S, I);
    if psnr ~= - 1
        D = imnoise(D, 'gaussian', 0, 10 ^ (- psnr / 10) * mean(mean(D .^ 2)));
    end
    K = ProtoPowerExpansionCoefficientsDB(n + ceil(n / 2) + 1, - length(D) - 50, length(D) + 50);
    M = zeros(2 * (ceil(n / 2) + 1) ^ 2, nchoosek(n + 2, 2));
    row = 1;
    for r = 0 : ceil(n / 2)
        for s = 0 : ceil(n / 2)
            col = 1;
            for i = 0 : n
                for j = 0 : n - i
                    if 0 <= i + r - 1
                        M(row + 0, col) = (i + r) * sum(sum(K(j + s + 1, t : t + length(D) - 1)' * K(i + r + 0, t : t + length(D) - 1) .* D));
                    else
                        M(row + 0, col) = 0;
                    end
                    if 0 <= j + s - 1
                        M(row + 1, col) = (j + s) * sum(sum(K(j + s + 0, t : t + length(D) - 1)' * K(i + r + 1, t : t + length(D) - 1) .* D));
                    else
                        M(row + 1, col) = 0;
                    end
                    col = col + 1;
                end
            end
            row = row + 2;
        end
    end
    C = quadprog(M' * M, [], [], [], eye(1, size(M, 2)), 1, [], [], [], optimset('display', 'off'))';
end