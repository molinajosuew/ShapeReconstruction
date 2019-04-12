






%%%%%%%%%%%%%%
% pixelizing and all-white image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('      --- Generating the pixels of an all-white image ----')
tic

% Applying the spline PSF on the binary image, followed by sampling
white_pixs      = pixelization(ones(length(y_range) , length(x_range)) , x_range , y_range , spline_degree);

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')




%%%%%%%%%%%%%%%%%%%%%%%%%
% finding the moments of an all-white image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('      --- Generalized moments of an all-white image ----')
tic

% Estimating the generalized moments of the image
[aux1 , aux2 , white_moments]   = moment_generator(white_pixs , m_x , m_y , coef_g , coef_dg);

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forming the inequalities regarding the polynomial coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('      --- Forming the inequality matrix ----')
tic

% Forming the matrix storing the inequality coefficients
even_moment_matrix  = matrix4inequality(sum_degree , moments_x_y , white_moments);

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')








%%%%%%%%%%%
% quadratic programming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('      --- Solving the quadratic programming ----')
tic

% initializing the minimization with a centered circle
x0              = [36 0 0 -1 0 -1 0 0 0 0 0 0 0 0 0]';

% minimization settings
options         = optimset('MaxIter' , 1e8);

% minimization problem with inequality and equality constraints
poly_coefs_QP   = quadprog(moment_matrix' * moment_matrix  ,  zeros((sum_degree+1)*(sum_degree+2)/2,1)  ,  ...
                           even_moment_matrix  ,  zeros(size(even_moment_matrix,1),1)  ,  ...
                           [1 , zeros(1,(sum_degree+1)*(sum_degree+2)/2-1)]  ,  [1] , ...
                           [] , [] , x0 , options);                    
                       
ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')




%%%%%%%%%%%%%%%%
% forming the estimated curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('      --- Forming the estimated AC image ----')
tic

binary_image_QP = AC_image(poly_coefs_QP , x_range , y_range);

figure,
imshow(binary_image_QP)

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')



%%%%%%%%%%
% Saving the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('      --- Saving the estimated AC image ----')
tic

imwrite(binary_image_QP , ['./Results/QP_BinaryImage_SNR' , num2str(Pixel_SNR) , '.tif'])

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')