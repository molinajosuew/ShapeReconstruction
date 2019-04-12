% This subroutine aims at generating the binary image of an algebraic curve
% and converting it into a set of pixels

disp('      --- Generating the AC image ----')
tic

% Generating the binary image of the algebraic curve
binary_image    = AC_image(poly_coefs , x_range , y_range);

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')


%%%%%%%%%%%%%%


disp('      --- Generating the pixels ----')
tic

% Applying the spline PSF on the binary image, followed by sampling
[pixelized , m_x , m_y]     = pixelization(binary_image , x_range , y_range , spline_degree);

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')


%%%%%%%%%%%%%%


disp('      --- Saving the results ----')
tic

% Saving the results
imwrite(binary_image , './Results/Clean_BinaryImage.tif')
imwrite(imresize(pixelized , round(length(binary_image) / length(pixelized)) , 'nearest') , ...
        './Results/Clean_Pixelized.tif')

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')