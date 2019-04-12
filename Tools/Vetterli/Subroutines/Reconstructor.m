% This subroutine aims at estimating the coefficients of the polynomial

disp('      --- Evaluating the generalized moments ----')
tic

% Estimating the generalized moments of the image
[moments_dx_y , moments_x_dy , moments_x_y]   = moment_generator(pixelized_noisy , m_x , m_y , coef_g , coef_dg);

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')


%%%%%%%%%%%%%%


disp('      --- Forming the null-producing matrix ----')
tic

% Forming the matrix, the null-space of which contains the coefficients
moment_matrix   = matrix4nullspace(sum_degree , moments_dx_y , moments_x_dy);

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')


%%%%%%%%%%%%%%


% disp('      --- Cadzow denoising ----')
% tic
% 
% % denoising the matrix
% moment_matrix   = Cadzow(moment_matrix , 1 , 100);
% 
% ttt             = toc;
% disp(['          It took ' , num2str(ttt) , ' seconds'])
% disp(' ')


%%%%%%%%%%%%%%


disp('      --- Row normalization ----')
tic

% normalizing the rows of the moment_matrix
for row_ind = 1 : size(moment_matrix , 1)
    moment_matrix(row_ind , :)  = moment_matrix(row_ind , :) / norm(moment_matrix(row_ind , :));
end

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')


%%%%%%%%%%%%%%


disp('      --- Finding the best match for the null-vector ----')
tic

% finding the null-space
[U,S,V]         = svd(moment_matrix);
poly_coefs_est  = V(: , end)';
binary_image_est= AC_image(poly_coefs_est , x_range , y_range);

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')


%%%%%%%%%%%%%%


disp('      --- Forcing partial consistency ----')
tic

% resizing the pixels
BI_est_small    = imresize(binary_image_est , size(pixelized_noisy));

% checking whether black and white should be reversed
if sum(sum( (BI_est_small - 0.5) .* pixelized_noisy )) < 0
    poly_coefs_est  = -poly_coefs_est;
    binary_image_est= 1 - binary_image_est;
end

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')


%%%%%%%%%%%%%%


disp('      --- Saving the result ----')
tic

% Saving the results
imwrite(binary_image_est , ['./Results/Est_BinaryImage_SNR' , num2str(Pixel_SNR) , '.tif'])

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')


%%%%%%%%%%%%%%


disp('      --- Showing the differences ----')
tic

% plotting the difference
figure,
imshow(abs(binary_image_est - binary_image))
title('Difference between the reconstructed (null-space method) and original images')
xlabel(['PSNR = ' , num2str( -10*log10( mean(mean((binary_image_est - binary_image).^2)) ) )])

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')



%%%%%%%%%%%%%%


if 1
    disp('      --- Consistency with the sign of the pixels ----')
    tic
    
    % Taking into account the sign of the pixels
    Sign_Consistency
    
    ttt             = toc;
    disp(['          It took ' , num2str(ttt) , ' seconds'])
    disp(' ')



    %%%%%%%%%%%%%%



    disp('      --- Postprocessing for Enhancement ----')
    tic
    
    % enhancing the outcome by quadratic programming and gradient movement
    Enhancement
    
    ttt             = toc;
    disp(['          It took ' , num2str(ttt) , ' seconds'])
    disp(' ')
end