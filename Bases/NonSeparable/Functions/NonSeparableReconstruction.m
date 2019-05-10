function [C, D] = NonSeparableReconstruction(I, n, L, psnr)
    m_x = (size(I, 2) - 1) / L;
    m_y = (size(I, 1) - 1) / L;

    if floor(m_y) ~= m_y || floor(m_x) ~= m_x
        disp("image domain and resolution inconsistency");
        return;
    end
    
    D_t = conv2(wjmBSpline(2 * n - 1, m_y), wjmBSpline(2 * n - 1, m_x), I) / (m_x * m_y);
    D_t_x = (size(D_t, 2) - size(I, 2)) / 2 + 1;
    D_t_y = (size(D_t, 1) - size(I, 1)) / 2 + 1;

    D = D_t(D_t_y - floor((D_t_y - 1) / m_y) * m_y : m_y : end, D_t_x - floor((D_t_x - 1) / m_x) * m_x : m_x : end);
    D_x = (size(D, 2) - size(I(1 : m_y : end, 1 : m_x : end), 2)) / 2 + 1;
    D_y = (size(D, 1) - size(I(1 : m_y : end, 1 : m_x : end), 1)) / 2 + 1;
    
    if psnr ~= - 1
        D = imnoise(D, 'gaussian', 0, 10 ^ (- psnr / 10) * mean(mean(D .^ 2)));
    end
    
    K = NonSeparableExpansionCoefficients(2 * n - 1, 2 * L);
    K_c = (size(K, 4) + 1) / 2;
    K_x = K_c + 1 - D_x : K_c + size(D, 2) - D_x;
    K_y = K_c + 1 - D_y : K_c + size(D, 1) - D_y;

    A = zeros(2 * (ceil(n / 2) + 1) ^ 2, nchoosek(n + 2, 2));
    row = 1;

    for a = 0 : ceil(n / 2)
        for b = 0 : ceil(n / 2)
            col = 1;

            for i = 0 : n
                for j = 0 : n - i
                    A(row + 0, col) = (NonSeparableMoment(n + a + b - 1, i + a - 1, j + b + 0, D, K, K_x, K_y) - NonSeparableMoment(n + a + b - 1, i + a, j + b, D, K, K_x, K_y)) * n * NonSeparableFactor(n, a, b, i, j, 2 * L) / (2 * L);
                    A(row + 1, col) = (NonSeparableMoment(n + a + b - 1, i + a + 0, j + b - 1, D, K, K_x, K_y) - NonSeparableMoment(n + a + b - 1, i + a, j + b, D, K, K_x, K_y)) * n * NonSeparableFactor(n, a, b, i, j, 2 * L) / (2 * L);

                    col = col + 1;
                end
            end

            row = row + 2;
        end
    end
    
    C = quadprog(100000 * (A' * A), [], [], [], eye(1, size(A, 2)), 1, [], [], [], optimset('MaxIter', 200, 'Algorithm', 'interior-point-convex', 'display', 'off'))';
end