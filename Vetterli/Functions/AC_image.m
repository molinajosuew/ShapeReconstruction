function binary_image = AC_image(poly_coefs , x_range , y_range)
%
% binary_image = AC_image(poly_coefs , x_range , y_range)
%
% This function generates the binary image associated with the 
% characteristic image of an algebraic curve (AC).
%
% "poly_coefs" is a row vector that contains the coefficients of the
%              two-variable polynomial representing the algebraic curve. If
%                               p(x,y) = 1 + 2y - x + 3xy
%              denotes a sample polynomial, the ordering of "poly_coefs" is
%              such that:
%                           "poly_coefs" = [1 , 2 , -1 , 3].
%
% "x_range" and "y_range" are row vectors that reveal the grid for
%           generating the images.
%
%
% "binary_image" is the 2D binary matrix that represents the characteristic
%                image of the algebraic curve. Its number of rows and 
%                columns are equal to the lengths of "y_range" and 
%                "x_range", respectively.



% finding the polynomial degree
sum_degree  = ceil(sqrt(2 * length(poly_coefs) + 1/4) - 3/2);


% finding all x^i and y^j
xx_allPower = ones(sum_degree+1 , length(x_range));
yy_allPower = ones(sum_degree+1 , length(y_range));
for pp = 1 : sum_degree
    xx_allPower(pp+1 , :)   = xx_allPower(pp , :) .* x_range;
    yy_allPower(pp+1 , :)   = yy_allPower(pp , :) .* y_range(end:-1:1);
end




% auxillary variables to accomodate for the polynomial implementation
x_power     = 0;
y_power     = 0;


% Evaluating the polynomial over the given range
binary_image= zeros(length(y_range) , length(x_range));

term        = 0;
while term < length(poly_coefs)
    term    = term + 1;
    
    binary_image    = binary_image + poly_coefs(term) * ...
                                    (yy_allPower(y_power+1 , :)' * xx_allPower(x_power+1 , :));
    
    if y_power == 0
        y_power     = x_power + y_power + 1;
        x_power     = 0;
        
    else
        y_power     = y_power - 1;
        x_power     = x_power + 1;
    end

end



% Defining the binary image
binary_image(binary_image <= 0) = 0;
binary_image(binary_image >  0) = 1;



% plotting the binary image
if 0
    figure,
    imshow(binary_image)
end