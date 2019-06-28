function [QP_E, QP, pixelized_noisy] = VetterliReconstruction(binary_image, x_a, x_b, y_a, y_b, x_n, y_n, psnr)
    x_range = linspace(x_a, x_b, x_n);
    y_range = linspace(y_a, y_b, y_n);
    [pixelized, m_x, m_y] = pixelization(binary_image, x_range, y_range, 4);

    if psnr ~= - 1
        pixelized_noisy = imnoise(pixelized, 'gaussian', 0, 10 ^ (- psnr / 10) * mean(mean(pixelized .^ 2)));
    else
        pixelized_noisy = pixelized;
    end

    load(['coefs_n_4_m_7_maxPower_', num2str(6), '_xMax_16.mat']);
    [moments_dx_y, moments_x_dy, ~] = moment_generator(pixelized_noisy, m_x, m_y, coef_g, coef_dg);
    moment_matrix = matrix4nullspace(4, moments_dx_y, moments_x_dy);

    for row_ind = 1 : size(moment_matrix, 1)
        moment_matrix(row_ind, :) = moment_matrix(row_ind, :) / norm(moment_matrix(row_ind, :));
    end

    x_power = 0;
    y_power = 0;
    x_power_vec = zeros(1, 15);
    y_power_vec = zeros(1, 15);
    term = 0;

    while x_power + y_power <= 4
        term = term + 1;
        x_power_vec(term) = x_power;
        y_power_vec(term) = y_power;

        if y_power == 0
            y_power = x_power + y_power + 1;
            x_power = 0;
        else
            y_power = y_power - 1;
            x_power = x_power + 1;
        end
    end

    Poly_eval_mat = zeros(length(m_x) * length(m_y), 15);

    for m_x_ind = 1 : length(m_x)
        for m_y_ind = 1 : length(m_y)
            row_ind = m_y_ind + (m_x_ind - 1) * length(m_y);
            Poly_eval_mat(row_ind, :) = (m_x(m_x_ind) .^ x_power_vec) .* (m_y(end - m_y_ind + 1) .^ y_power_vec);
        end
    end

    HH = moment_matrix' * moment_matrix;
    ff = zeros(15, 1);
    A_eq = [1, zeros(1, 14)];
    b_eq = sign(pixelized_noisy(m_y == 0, m_x == 0) - 0.5);
    A_ineq = Poly_eval_mat;
    b_ineq = ones(length(m_x) * length(m_y), 1) * inf;

    for m_x_ind = 1 : length(m_x)
        for m_y_ind = 1 : length(m_y)
            if pixelized_noisy(m_y_ind, m_x_ind) >= min(0.7, max(max(pixelized_noisy)))
                row_ind = m_y_ind + (m_x_ind - 1) * length(m_y);
                b_ineq(row_ind) = 0;
                A_ineq(row_ind, :) = - A_ineq(row_ind, :);
            elseif pixelized_noisy(m_y_ind, m_x_ind) <= 0.3
                row_ind = m_y_ind + (m_x_ind - 1) * length(m_y);
                b_ineq(row_ind) = 0;
            end
        end
    end

    options = optimset('MaxIter', 200, 'Algorithm', 'interior-point-convex', 'display', 'off');
    poly_coefs_QP = quadprog(10000 * HH, ff, A_ineq(b_ineq ~= inf, :), b_ineq(b_ineq ~= inf), A_eq, b_eq, [], [], [], options);
    binary_image_QP = AC_image(poly_coefs_QP, x_range, x_range);
    QP = binary_image_QP;
    pixelized_QP = pixelization(binary_image_QP , x_range , y_range , 4);
    maxIter = 20;
    delta = 1e-4;
    current_Err = 1;
    new_Err = 0;
    iter = 0;

    while (((current_Err - new_Err) / current_Err) > 0.001 ) && (iter < maxIter)
        iter = iter + 1;
        gradient_matrix = zeros(length(m_x) * length(m_y), 15);

        for coef_ind = 1 : 15
            new_coef = poly_coefs_QP;
            new_coef(coef_ind) = new_coef(coef_ind) + delta;
            binary_image_new = AC_image(new_coef , x_range, y_range);
            pixelized_new = pixelization(binary_image_new, x_range, y_range, 4);
            aux = (pixelized_new - pixelized_QP) / delta;
            gradient_matrix(:, coef_ind) = aux(:);
        end

        diff = pixelized_QP - pixelized_noisy;
        diff = diff(:);
        current_Err = norm(diff);
        change = - pinv(gradient_matrix) * diff;
        binary_image_QP = AC_image(poly_coefs_QP + change, x_range, y_range);
        pixelized_QP = pixelization(binary_image_QP, x_range, y_range, 4);
        new_Err = sum(sum((pixelized_QP - pixelized_noisy) .^ 2)) .^ 0.5;

        while new_Err > current_Err * 1.001
            change = change / 2;
            binary_image_QP = AC_image(poly_coefs_QP + change, x_range, y_range);
            pixelized_QP = pixelization(binary_image_QP, x_range, y_range, 4);
            new_Err = sum(sum((pixelized_QP - pixelized_noisy) .^ 2)) .^ 0.5;
        end

        poly_coefs_QP = poly_coefs_QP + change;
    end

    QP_E = binary_image_QP;
end