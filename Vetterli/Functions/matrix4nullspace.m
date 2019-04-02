function moment_matrix = matrix4nullspace(sum_degree , moments_dx_y , moments_x_dy)
%
% moment_matrix = matrix4nullspace(sum_degree , moments_dx_y , moments_x_dy)
%
% We assume "moments_dx_y" and "moments_x_dy" are the moments of the
% indicator image of an algebraic curve (AC) with maximum degree sum of
% "sum_degree". It can be shown that the polynomial coefficients of the
% algebraic curve responsible for the moments, falls in the null space of a
% matrix formed by the moments. The aim of this function is to form this
% matrix.
%
% "sum_degree" is the maximum degree of the 2-variate polynomial as a sum.
%              In other words, if 
%                         P(x,y) = \sum_{i,j >= 0} c(i,j) x^i y^j
%              denotes the 2-variate polynomial, "sum_degree" shall be the
%              maximum of 'i+j' over all paird of (i,j) such that c(i,j) is
%              non-zero.
%
% "moments_dx_y" is a matrix, the (i,j) element of which contains the
%                inner product of the original AC indicator image with
%                       d/dx (x^(i-1) * y^(j-1) * h(x) * g(y)) 
%                where 'h(.)' and 'g(.)' are some given functions.
%
% "moments_x_dy" is a matrix, the (i,j) element of which contains the
%                inner product of the original AC indicator image with
%                       d/dy (x^(i-1) * y^(j-1) * g(x) * h(y))
%                where 'h(.)' and 'g(.)' are some given functions.
%
% "moment_matrix" is the aforementioned matrix of moments in a suitable
%                 order such that it contains the vector of polynomial
%                 coefficients in its null-space.



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
% Generating the first part of the matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% auxillary vafriables for determining the number of rows
xx_max              = size(moments_dx_y , 1) - sum_degree - 1;
yy_max              = size(moments_dx_y , 2) - sum_degree - 1;

if xx_max < -1
    xx_max          = -1;
end

if yy_max < -1
    yy_max          = -1;
end

% defining the first part
moment_matrix_dx_y  = zeros( (xx_max+1)*(yy_max+1)  ,  (sum_degree+1) * (sum_degree+2) / 2 );
row_ind             = 0;
for xx_ind = 0 : xx_max
    for yy_ind = 0 : yy_max
        row_ind     = row_ind + 1;
        for col_ind = 1 : (sum_degree+1) * (sum_degree+2) / 2
            moment_matrix_dx_y(row_ind , col_ind)   = moments_dx_y(x_power_vec(col_ind) + xx_ind + 1  ,  y_power_vec(col_ind) + yy_ind + 1);
        end
    end
end





%%%%%%%%%%%%%%%%%%%%
% Generating the second part of the matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% auxillary vafriables for determining the number of rows
xx_max              = size(moments_x_dy , 1) - sum_degree - 1;
yy_max              = size(moments_x_dy , 2) - sum_degree - 1;

if xx_max < -1
    xx_max          = -1;
end

if yy_max < -1
    yy_max          = -1;
end

% defining the second part
moment_matrix_x_dy  = zeros( (xx_max+1)*(yy_max+1)  ,  (sum_degree+1) * (sum_degree+2) / 2 );
row_ind             = 0;
for xx_ind = 0 : xx_max
    for yy_ind = 0 : yy_max
        row_ind     = row_ind + 1;
        for col_ind = 1 : (sum_degree+1) * (sum_degree+2) / 2
            moment_matrix_x_dy(row_ind , col_ind)   = moments_x_dy(x_power_vec(col_ind) + xx_ind + 1  ,  y_power_vec(col_ind) + yy_ind + 1);
        end
    end
end






%%%%%%%%%%%%%%
% combining the two matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
moment_matrix       = [moment_matrix_dx_y ; moment_matrix_x_dy];

% normalizing the rows
if 0
    for row_ind = 1 : size(moment_matrix , 1)
        moment_matrix(row_ind , :)  = moment_matrix(row_ind , :) / norm(moment_matrix(row_ind , :));
    end
end