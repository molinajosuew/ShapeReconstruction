clear;
clc;
x_n = 512;
y_n = 512;
x_range = linspace(- 5, 5, x_n);
y_range = linspace(- 5, 5, y_n);
binary_image = GetImageOfPowerPolynomial(GetRandomPowerPolynomial(4, - 5, 5, - 5, 5, false), - 5, 5, - 5, 5, x_n, y_n);
% binary_image = RandomStablyBounded4thDegreeBivariatePolynomial(- 5, 5, - 5, 5, x_n, y_n);
psnr = - 1;
[pixelized, m_x, m_y] = pixelization(binary_image, x_range, y_range, 4);
if psnr ~= - 1
    pixelized_noisy = imnoise(pixelized, 'gaussian', 0, 10 ^ (- psnr / 10) * mean(mean(pixelized .^ 2)));
else
    pixelized_noisy = pixelized;
end
load(['coefs_n_4_m_7_maxPower_', num2str(6), '_xMax_16.mat']);
[moments_dx_y, moments_x_dy, ~] = moment_generator(pixelized_noisy, m_x, m_y, coef_g, coef_dg);
moment_matrix = matrix4nullspace(4, moments_dx_y, moments_x_dy);
[~, ~, V] = svd(moment_matrix);
poly_coefs_est = V(:, end)';
binary_image_est = AC_image(poly_coefs_est, x_range, x_range);
BI_est_small = imresize(binary_image_est, size(pixelized_noisy));
if sum(sum((BI_est_small - 0.5) .* pixelized_noisy)) < 0
    poly_coefs_est  = - poly_coefs_est;
    binary_image_est= 1 - binary_image_est;
end
x_power_vec = zeros(1, 15);
y_power_vec = zeros(1, 15);
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
b_eq = sign( pixelized_noisy(m_y == 0, m_x == 0) - 0.5);
x0 = poly_coefs_est';
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
options = optimset('MaxIter', 2e2, 'Algorithm', 'interior-point-convex');
[poly_coefs_QP, cost, exitflag] = quadprog(HH * 10000, ff, A_ineq(b_ineq ~= inf, :), b_ineq(b_ineq ~= inf), A_eq, b_eq, [], [], x0, options);
exitflag;
binary_image_QP = AC_image(poly_coefs_QP, x_range, x_range);
pixelized_QP    = pixelization(binary_image_QP , x_range , y_range , 4);
% imshow(abs(binary_image - binary_image_QP));
maxIter         = 20;
delta           = 1e-4;
current_Err     = 1;
new_Err         = 0;
iter            = 0;
while ( ((current_Err - new_Err) / current_Err) > 0.001 ) && ( iter < maxIter )
    iter        = iter+1;
    gradient_matrix         = zeros(length(m_x)*length(m_y)  ,  15);
    for coef_ind = 1 : 15
        new_coef            = poly_coefs_QP;
        new_coef(coef_ind)  = new_coef(coef_ind) + delta;
        binary_image_new    = AC_image(new_coef , x_range , y_range);
        pixelized_new       = pixelization(binary_image_new , x_range , y_range , 4);
        aux                 = (pixelized_new - pixelized_QP) / delta;
        gradient_matrix(: , coef_ind)   = aux(:);
    end
    diff    = pixelized_QP - pixelized_noisy;
    diff    = diff(:);
    current_Err     = norm(diff);
    change  = - pinv(gradient_matrix) * diff;
    binary_image_QP = AC_image(poly_coefs_QP + change , x_range , y_range);
    pixelized_QP    = pixelization(binary_image_QP , x_range , y_range , 4);
    new_Err = sum(sum( (pixelized_QP-pixelized_noisy).^2 )).^0.5;
    while new_Err > current_Err * 1.001
        change  = change / 2;
        binary_image_QP = AC_image(poly_coefs_QP + change , x_range , y_range);
        pixelized_QP    = pixelization(binary_image_QP , x_range , y_range , 4);
        new_Err = sum(sum( (pixelized_QP-pixelized_noisy).^2 )).^0.5;
    end
    poly_coefs_QP   = poly_coefs_QP + change;
    if 1
        P_SNR       = 10 * log10( mean(mean(pixelized_noisy.^2)) / mean(mean( (pixelized_QP-pixelized_noisy).^2 )));
        disp(['      >>> pixel consistency: ' , num2str(P_SNR) , ' dB , maximum error = ' , num2str(max(max( abs(pixelized_QP-pixelized_noisy) )))])
    end
end
imshow(abs(binary_image_QP - binary_image))