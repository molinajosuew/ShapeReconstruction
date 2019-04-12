function [R, J2] = ProtoPowerReconstruction1(I, n, x_a, x_b, y_a, y_b)
    m_x = (size(I, 2) - 1) / (x_b - x_a);
    m_y = (size(I, 1) - 1) / (y_b - y_a);

    if floor(m_y) ~= m_y || floor(m_x) ~= m_x
        disp("image domain and resolution inconsistency"); % This can be fixed by up-sampling and shifting splines, I think.
        return;
    end
    
    J1 = conv2(wjmBSpline(n + ceil(n / 2), m_y), wjmBSpline(n + ceil(n / 2), m_x), I) / (m_x * m_y);
    J1_rx = (size(J1, 2) - size(I, 2)) / 2 + 1;
    J1_ry = (size(J1, 1) - size(I, 1)) / 2 + 1;

    J2 = J1(J1_ry - floor((J1_ry - 1) / m_y) * m_y : m_y : end, J1_rx - floor((J1_rx - 1) / m_x) * m_x : m_x : end);
    J2_rx = (size(J2, 2) - size(I(1 : m_y : end, 1 : m_x : end), 2)) / 2 + 1;
    J2_ry = (size(J2, 1) - size(I(1 : m_y : end, 1 : m_x : end), 1)) / 2 + 1;

    l = max(abs([1 - J2_rx + x_a, 1 - J2_ry + y_a, size(J2, 2) - J2_rx + x_a, size(J2, 1) - J2_ry + y_a]));
    
    C = ProtoPowerCoefficients(n + ceil(n / 2), - l, l);
    C_r = (size(C, 2) + 1) / 2;
    C_x = C_r + 1 - J2_rx + x_a : C_r + size(J2, 2) - J2_rx + x_a;
    C_y = C_r + 1 - J2_ry + y_a : C_r + size(J2, 1) - J2_ry + y_a;

    M = zeros(2 * (ceil(n / 2) + 1) ^ 2, nchoosek(n + 2, 2));
    row = 1;

    for r = 0 : ceil(n / 2)
        for s = 0 : ceil(n / 2)
            col = 1;

            for i = 0 : n
                for j = 0 : n - i
                    if 0 <= i + r - 1 && 0 <= j + s + 0
                        M(row + 0, col) = (i + r) * sum(sum(C(j + s + 1, C_y)' * C(i + r + 0, C_x) .* J2));
                    else
                        M(row + 0, col) = 0;
                    end

                    if 0 <= i + r + 0 && 0 <= j + s - 1
                        M(row + 1, col) = (j + s) * sum(sum(C(j + s + 0, C_y)' * C(i + r + 1, C_x) .* J2));
                    else
                        M(row + 1, col) = 0;
                    end

                    col = col + 1;
                end
            end

            row = row + 2;
        end
    end

%     R = quadprog(M' * M, [], [], [], eye(1, size(M, 2)), 1)';
    R = lsqlin(M, zeros(size(M, 1), 1), [], [], eye(1, size(M, 2)), 1, [], [])';
end