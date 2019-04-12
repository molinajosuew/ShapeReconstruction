% This subroutine emforces sign consistency with the available pixels; in
% simple words, it tries to find the closest fit to the null space of the
% moments that generates almost the same pixels





%%%%%%%%%%%%%%%%%%%%%
% Finding the polynomial degree vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% auxillary variables for storing the monomial powers
x_power     = 0;
y_power     = 0;

x_power_vec = zeros(1 , (sum_degree+1) * (sum_degree+2) / 2 );
y_power_vec = zeros(1 , (sum_degree+1) * (sum_degree+2) / 2 );



% finding the order of the monomials up to degree "sum_degree"
term        = 0;
while x_power + y_power <= sum_degree
    term    = term + 1;
    
    x_power_vec(term)   = x_power;
    y_power_vec(term)   = y_power;
    
    if y_power == 0
        y_power     = x_power + y_power + 1;
        x_power     = 0;
        
    else
        y_power     = y_power - 1;
        x_power     = x_power + 1;
        
    end

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polynomial evaluation matrix over the sampling grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Poly_eval_mat   = zeros( length(m_x)*length(m_y)  ,  (sum_degree+1) * (sum_degree+2) / 2 );
for m_x_ind = 1 : length(m_x)
    for m_y_ind = 1 : length(m_y)
        row_ind     = m_y_ind + (m_x_ind - 1) * length(m_y);
        Poly_eval_mat(row_ind , :)  = (m_x(m_x_ind) .^ x_power_vec) .* (m_y(end-m_y_ind+1) .^ y_power_vec);
    end
end





%%%%%%%%%%%%%%%%%%%%%%
% Quadratic programming for initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HH              = moment_matrix' * moment_matrix;
ff              = zeros((sum_degree+1) * (sum_degree+2) / 2   ,   1);

A_eq            = [1   ,   zeros(1 , (sum_degree+1) * (sum_degree+2) / 2 - 1)];
b_eq            = sign( pixelized_noisy(m_y==0 , m_x==0) - 0.5 );

x0              = poly_coefs_est';

A_ineq          = Poly_eval_mat;
b_ineq          = ones(length(m_x)*length(m_y) , 1) * inf;
for m_x_ind = 1 : length(m_x)
    for m_y_ind = 1 : length(m_y)
        if pixelized_noisy(m_y_ind , m_x_ind) >= min(0.7 , max(max(pixelized_noisy)))
            row_ind     = m_y_ind + (m_x_ind - 1) * length(m_y);
            b_ineq(row_ind)         = 0;
            A_ineq(row_ind , :)     = -A_ineq(row_ind , :);
            
        elseif pixelized_noisy(m_y_ind , m_x_ind) <= 0.3
            row_ind     = m_y_ind + (m_x_ind - 1) * length(m_y);
            b_ineq(row_ind)         = 0;
            
        end
    end
end

options         = optimset('MaxIter' , 2e2 , 'Algorithm' , 'interior-point-convex');


% The initial minimization, the output of which is used as the initial
% point for the next iteration
[poly_coefs_QP , cost , exitflag]   = quadprog(HH*10000 , ff , A_ineq(b_ineq~=inf , :) , b_ineq(b_ineq~=inf) , A_eq , b_eq , [] , [] , x0 , options);
%[poly_coefs_QP , cost , exitflag]   = quadprog(HH , ff , A_ineq , b_ineq , A_eq , b_eq , [] , [] , [] , options);
exitflag;


% finding the binary image of the current AC estimate
binary_image_QP = AC_image(poly_coefs_QP , x_range , y_range);

% computing the associated pixels
pixelized_QP    = pixelization(binary_image_QP , x_range , y_range , spline_degree);

if 1
    P_SNR       = 10 * log10( mean(mean(pixelized_noisy.^2)) / mean(mean( (pixelized_QP-pixelized_noisy).^2 )));
    disp(['      >>> pixel consistency: ' , num2str(P_SNR) , ' dB , maximum error = ' , num2str(max(max( abs(pixelized_QP-pixelized_noisy) )))])
end






%%%%%%%%%%%%%%
% Checking the QP result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotting the difference of the output with the original image
figure,
imshow(abs(binary_image_QP - binary_image))
title('Difference between the reconstructed (zero-Thresh QP) and original images')
xlabel(['PSNR = ' , num2str( -10*log10( mean(mean((binary_image_QP - binary_image).^2)) ) )])

drawnow
