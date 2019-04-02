




%%%%%%%%%%%%%%%%%%
% Iterations to improve consistency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxIter         = 20;
delta           = 1e-4;     % dx for approximating the gradient

current_Err     = 1;
new_Err         = 0;
iter            = 0;

% main iterations
while ( ((current_Err - new_Err) / current_Err) > 0.001 ) && ( iter < maxIter )

    iter        = iter+1;
    
    % the gradient of the coefficient with respect to the pixels
    gradient_matrix         = zeros(length(m_x)*length(m_y)  ,  (sum_degree+1) * (sum_degree+2) / 2);
    for coef_ind = 1 : (sum_degree+1) * (sum_degree+2) / 2
        %coef_ind
        
        % marginally changing the "coef_ind"th coefficient
        new_coef            = poly_coefs_QP;
        new_coef(coef_ind)  = new_coef(coef_ind) + delta;
    
        % finding the binary image for the modified coefficients
        binary_image_new    = AC_image(new_coef , x_range , y_range);

        % computing the associated pixels
        pixelized_new       = pixelization(binary_image_new , x_range , y_range , spline_degree);
    

        aux                 = (pixelized_new - pixelized_QP) / delta;
        gradient_matrix(: , coef_ind)   = aux(:);
        
    end
        
    
    % evaluating the pixelization consistency
    diff    = pixelized_QP - pixelized_noisy;
    diff    = diff(:);
    
    current_Err     = norm(diff);
   
    
    % possible modification in coefficients
    change  = - pinv(gradient_matrix) * diff;
    
    
    
    
    %%%%%%%%%%%%%
    % tuning the step size
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % finding the binary image of the new AC estimate
    binary_image_QP = AC_image(poly_coefs_QP + change , x_range , y_range);
    
    % computing the associated pixels
    pixelized_QP    = pixelization(binary_image_QP , x_range , y_range , spline_degree);
    
    % evaluating the new error
    new_Err = sum(sum( (pixelized_QP-pixelized_noisy).^2 )).^0.5;
    
    % decreasing the step size if the error has increased
    while new_Err > current_Err * 1.001
        change  = change / 2;
        
        % finding the binary image of the new AC estimate
        binary_image_QP = AC_image(poly_coefs_QP + change , x_range , y_range);
        
        % computing the associated pixels
        pixelized_QP    = pixelization(binary_image_QP , x_range , y_range , spline_degree);
        
        % evaluating the new error
        new_Err = sum(sum( (pixelized_QP-pixelized_noisy).^2 )).^0.5;
    end
    
    % updating the coefficients
    poly_coefs_QP   = poly_coefs_QP + change;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%
    
    
    
    if 1
        P_SNR       = 10 * log10( mean(mean(pixelized_noisy.^2)) / mean(mean( (pixelized_QP-pixelized_noisy).^2 )));
        disp(['      >>> pixel consistency: ' , num2str(P_SNR) , ' dB , maximum error = ' , num2str(max(max( abs(pixelized_QP-pixelized_noisy) )))])
    end
    
    
end









%%%%%%%%%%%%%%
% Checking the final result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotting the difference of the output with the original image
figure,
imshow(abs(binary_image_QP - binary_image))
title('Difference between the reconstructed (QP+gradient) and original images')
xlabel(['PSNR = ' , num2str( -10*log10( mean(mean((binary_image_QP - binary_image).^2)) ) )])





%%%%%%%%%%%%
% Saving the final result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imwrite(binary_image_QP , ['./Results/QP_Est_BinaryImage_SNR' , num2str(Pixel_SNR) , '.tif'])
