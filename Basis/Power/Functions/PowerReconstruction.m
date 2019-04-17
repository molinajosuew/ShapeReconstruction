function [C, D] = PowerReconstruction(I, n, x_a, x_b, y_a, y_b, z)
    m_x = (size(I, 2) - 1) / (x_b - x_a);
    m_y = (size(I, 1) - 1) / (y_b - y_a);

    if floor(m_y) ~= m_y || floor(m_x) ~= m_x
        disp("image domain and resolution inconsistency"); % This can be fixed by up-sampling and shifting splines, I think.
        return;
    end
    
    D_t = conv2(wjmBSpline(n + ceil(n / 2), m_y), wjmBSpline(n + ceil(n / 2), m_x), I) / (m_x * m_y);
    D_t_rx = (size(D_t, 2) - size(I, 2)) / 2 + 1;
    D_t_ry = (size(D_t, 1) - size(I, 1)) / 2 + 1;

    D = D_t(D_t_ry - floor((D_t_ry - 1) / m_y) * m_y : m_y : end, D_t_rx - floor((D_t_rx - 1) / m_x) * m_x : m_x : end);
    D_rx = (size(D, 2) - size(I(1 : m_y : end, 1 : m_x : end), 2)) / 2 + 1;
    D_ry = (size(D, 1) - size(I(1 : m_y : end, 1 : m_x : end), 1)) / 2 + 1;
    
    if z ~= - 1
        D = awgn(D, z);
    end

    l = max(abs([1 - D_rx + x_a, 1 - D_ry + y_a, size(D, 2) - D_rx + x_a, size(D, 1) - D_ry + y_a]));
    
    K = PowerExpansionCoefficients(n + ceil(n / 2), - l, l);
    K_r = (size(K, 2) + 1) / 2;
    K_x = K_r + 1 - D_rx + x_a : K_r + size(D, 2) - D_rx + x_a;
    K_y = K_r + 1 - D_ry + y_a : K_r + size(D, 1) - D_ry + y_a;

    M = zeros(2 * (ceil(n / 2) + 1) ^ 2, nchoosek(n + 2, 2));
    row = 1;

    for r = 0 : ceil(n / 2)
        for s = 0 : ceil(n / 2)
            col = 1;

            for i = 0 : n
                for j = 0 : n - i
                    if 0 <= i + r - 1
                        M(row + 0, col) = (i + r) * sum(sum(K(j + s + 1, K_y)' * K(i + r + 0, K_x) .* D));
                    else
                        M(row + 0, col) = 0;
                    end

                    if 0 <= j + s - 1
                        M(row + 1, col) = (j + s) * sum(sum(K(j + s + 0, K_y)' * K(i + r + 1, K_x) .* D));
                    else
                        M(row + 1, col) = 0;
                    end

                    col = col + 1;
                end
            end

            row = row + 2;
        end
    end

    C = quadprog(M' * M, [], [], [], eye(1, size(M, 2)), 1)';
%     C = lsqlin(M, zeros(size(M, 1), 1), [], [], eye(1, size(M, 2)), 1, [], [])';
end