function denoised = Cadzow(noisy_matrix , null_dimension , max_iter)
%
% denoised = Cadzow(noisy_matrix , null_dimension , max_iter)
% 
% This function denoises a given matrix that have repeated elements and
% shall be rank deficient. The technique is usually referred to as Cadzow
% denoising and accounts for iterating between element-wise averaging and
% making it rank deficient.
%
% "noisy_matrix" is the input rectangular matrix which will be subject to
%                denoiding. It should have repeated elements, otherwise,
%                this method just wastes time.
%
% "null_dimension" is the hypothetical  dimension of column-wise null-space
%                  of the matrix. In each iteration of Cadzow, the smallest
%                  "null_dimension" singular values of the matrix will be
%                  made zero.
%
% "max_iter" is the number of iterations applied by the two mentioned
%            projections.
%
% "denoised" is the resulting matrix with the same size as "noisy_matrix".



denoised        = noisy_matrix;
for iter = 1 : max_iter
    % projecting the matrix onto the rank deficient matrices
    [U,S,V]         = svd(denoised);
    S(:  ,  end-null_dimension+1 : end)     = 0;
    denoised        = U*S*V';
            
    for row_ind = 1 : size(noisy_matrix , 1)
        for col_ind = 1 : size(noisy_matrix , 2)
            % preserving the repetetive structure
            rep_elements    = ( abs(noisy_matrix - noisy_matrix(row_ind , col_ind)) < 1e-9);
            denoised(rep_elements)  = mean( denoised(rep_elements) );
            
        end
    end
end