function [C, D] = ProtoNonSeparableBernsteinReconstruction(I, n, L, psnr)
    m_x = (size(I, 2) - 1) / L;
    m_y = (size(I, 1) - 1) / L;

    if floor(m_y) ~= m_y || floor(m_x) ~= m_x
        disp("image domain and resolution inconsistency"); % This can be fixed by up-sampling and shifting splines, I think.
        return;
    end
    
    D_t = conv2(wjmBSpline(2 * n - 1, m_y), wjmBSpline(2 * n - 1, m_x), I) / (m_x * m_y);
    D_t_rx = (size(D_t, 2) - size(I, 2)) / 2 + 1;
    D_t_ry = (size(D_t, 1) - size(I, 1)) / 2 + 1;

    D = D_t(D_t_ry - floor((D_t_ry - 1) / m_y) * m_y : m_y : end, D_t_rx - floor((D_t_rx - 1) / m_x) * m_x : m_x : end);
    D_rx = (size(D, 2) - size(I(1 : m_y : end, 1 : m_x : end), 2)) / 2 + 1;
    D_ry = (size(D, 1) - size(I(1 : m_y : end, 1 : m_x : end), 1)) / 2 + 1;
    
    if psnr ~= - 1
        D = imnoise(D, 'gaussian', 0, 10 ^ (- psnr / 10) * mean(mean(D .^ 2)));
    end
    
    K = ProtoNonSeparableBernsteinExpansionCoefficients(2 * n - 1, 2 * L);
    K_r = (size(K, 4) + 1) / 2;
    K_x = K_r + 1 - D_rx : K_r + size(D, 2) - D_rx;
    K_y = K_r + 1 - D_ry : K_r + size(D, 1) - D_ry;

    M = zeros(2 * (ceil(n / 2) + 1) ^ 2, nchoosek(n + 2, 2));
    row = 1;

    for a = 0 : ceil(n / 2)
        for b = 0 : ceil(n / 2)
            col = 1;

            for i = 0 : n
                for j = 0 : n - i
                    M(row + 0, col) = (ProtoNonSeparableBernsteinMoment(n + a + b - 1, i + a - 1, j + b + 0, D, K, K_x, K_y) - ProtoNonSeparableBernsteinMoment(n + a + b - 1, i + a, j + b, D, K, K_x, K_y)) * n * ProtoNonSeparableBernsteinFactor(n, a, b, i, j, 2 * L) / L / 2;
                    M(row + 1, col) = (ProtoNonSeparableBernsteinMoment(n + a + b - 1, i + a + 0, j + b - 1, D, K, K_x, K_y) - ProtoNonSeparableBernsteinMoment(n + a + b - 1, i + a, j + b, D, K, K_x, K_y)) * n * ProtoNonSeparableBernsteinFactor(n, a, b, i, j, 2 * L) / L / 2;

                    col = col + 1;
                end
            end

            row = row + 2;
        end
    end

    C = quadprog(M' * M, [], [], [], eye(1, size(M, 2)), 1)';
%     C = lsqlin(M, zeros(size(M, 1), 1), [], [], eye(1, size(M, 2)), 1, [], [])';
end