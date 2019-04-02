function even_moment_matrix = matrix4inequality(sum_degree , moments_x_y , white_moments)
%
% even_moment_matrix = matrix4inequality(sum_degree , moments_x_y , white_moments)
%
% We assume "moments_x_y" is the matrix of the moments of of the indicator 
% image of an algebraic curve (AC) with maximum degree sum of "sum_degree".
% It can be shown that the even-order moments imply certain inequalities
% between linear combination of the polynomial coefficients. The aim of 
% this function is to form the matrix that contains this linear
% combinations.
%
% "sum_degree" is the maximum degree of the 2-variate polynomial as a sum.
%              In other words, if 
%                         P(x,y) = \sum_{i,j >= 0} c(i,j) x^i y^j
%              denotes the 2-variate polynomial, "sum_degree" shall be the
%              maximum of 'i+j' over all paird of (i,j) such that c(i,j) is
%              non-zero.
%
% "moments_x_y" is a square matrix, the (i,j) element of which contains the
%               inner product of the original AC indicator image with
%                             x^(i-1) * y^(j-1) * g(x) * g(y)
%                where 'g(.)' is some non-negative function.
%
% "white_moments" is the same as "moments_x_y", however, instead of the 
%                 image an algebraic curve and all-white image of the same 
%                 size is used. In particular, the size of "white_moments"
%                 should match that of "moments_x_y".
%
% "even_moment_matrix" is the aforementioned matrix of even-order moments 
%                      in a suitable order such that if C representes the
%                      vector of polynomial coefficients, we shall have
%                               even_moment_matrix * C  <  0



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





%%%%%%%%%%%%%%%%%%%%
% Generating partial matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% auxillary variables for determining the number of rows
xy_max              = ceil( (length(moments_x_y) - sum_degree - 1) / 2);

if xy_max < -1
    xy_max          = -1;
end

% defining the moment-ordered matrix both for the curve and san all-whiteimage
curve_matrix        = zeros( (xy_max+1)^2  ,  (sum_degree+1) * (sum_degree+2) / 2 );
white_matrix        = zeros( (xy_max+1)^2  ,  (sum_degree+1) * (sum_degree+2) / 2 );

row_ind             = 0;
for xx_ind = 0 : xy_max
    for yy_ind = 0 : xy_max
        row_ind     = row_ind + 1;
        for col_ind = 1 : (sum_degree+1) * (sum_degree+2) / 2;
            curve_matrix(row_ind , col_ind)     = moments_x_y(  x_power_vec(col_ind) + 2*xx_ind + 1  ,  y_power_vec(col_ind) + 2*yy_ind + 1);
            white_matrix(row_ind , col_ind)     = white_moments(x_power_vec(col_ind) + 2*xx_ind + 1  ,  y_power_vec(col_ind) + 2*yy_ind + 1);
        end
    end
end







%%%%%%%%%%%%%%
% Forming the final matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
even_moment_matrix  = [-curve_matrix  ;  white_matrix - curve_matrix];

% normalizing the rows
if 0
    for row_ind = 1 : size(moment_matrix , 1)
        even_moment_matrix(row_ind , :) = even_moment_matrix(row_ind , :) / norm(even_moment_matrix(row_ind , :));
    end
end