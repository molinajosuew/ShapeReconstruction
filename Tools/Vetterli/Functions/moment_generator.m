
function [moments_dx_y , moments_x_dy , moments_x_y] = moment_generator(pixelized , m_x , m_y , coef_g , coef_dh)

%
% [moments_dx_y , moments_x_dy , moments_x_y] = moment_generator(pixelized , m_x , m_y , coef_g , coef_dh)
%
% It is assumed that the values stored within the matrix "pixelized" come
% from observing an image through a symmetric 2D kernel B(x,y); using the
% inner-product notation, this is equivalent to
%
%        pixelized(i,j)   =   < Image(x,y)  ,  B(x - m_x(i) , y - m_y(j)) >
%
% where "m_x" and "m_y" contain the shifts (grid) of the kernel. More
% specifically, "pixelized" is 2D matrix containing values between 0 and 1 
% that represent the pixels of the blurred image, subject to the PSF B(x,y).
%
% This function aims at finding the generalized moments of the image by
% linearly combining the pixels. The coefficients in this linear
% combination are provided in "coef_g" and "coef_dg". In particular, it is
% assumed that the functions 'g(.)' and 'h(.)' exist such that
%
%                    d/dy (x^i * y^j * g(x) * h(y)) = 
%
%            \sum_{m,n}   coef_g(i+2 , m) * coef_dh(j+2 , n) 
%                       * B(x - coef_g(1 , m)   ,   y - coef_dh(1 , n))
%
%
% Note that the first rows in "coef_g" and "coef_dh" are reserved for the 
% indices and not the coefficients behind B(.,.).
%
% 
% "moments_dx_y" is a matrix with 'size(coef_dh , 1)-1' rows and 
%                'size(coef_g , 1)-1' columns. The (i,j) element shall
%                contain the inner product of the original image with
%                       d/dx (x^(i-1) * y^(j-1) * h(x) * g(y)) 
%
% "moments_x_dy" is a matrix with 'size(coef_g , 1)-1' rows and 
%                'size(coef_dh , 1)-1' columns. The (i,j) element shall
%                contain the inner product of the original image with
%                       d/dy (x^(i-1) * y^(j-1) * g(x) * h(y))
%
% "moments_x_y" is a square matrix with size 'size(coef_g , 1)-1'. The 
%               (i,j) element shall contain the inner product of the 
%               original image with
%                             x^(i-1) * y^(j-1) * g(x) * g(y)
%               This matrix is not directly related to the story above,
%               however, plays a role in some restrcitive inequalities.



if (length(m_x) ~= size(pixelized , 2))||(length(m_y) ~= size(pixelized , 1))
    m_x         = -floor( size(pixelized , 2)/2 ) : size(pixelized , 2)-floor( size(pixelized , 2)/2 )-1;
    m_y         = -floor( size(pixelized , 1)/2 ) : size(pixelized , 1)-floor( size(pixelized , 1)/2 )-1;
end
    
    

%%%%%%%
% Generating "moments_dx_y"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% making the size of the matrices consistent (X direction)
aux1            = coef_dh(1 , :);
aux2            = min(aux1(1) , m_x(1)) : max(aux1(end) , m_x(end));
coef_dh_new     = [zeros(size(coef_dh , 1)   , aux1(1) - aux2(1))   ,    coef_dh   ,   zeros(size(coef_dh , 1)   , aux2(end) - aux1(end))];
pixelized_new   = [zeros(size(pixelized , 1) , m_x(1) - aux2(1) )   ,   pixelized  ,   zeros(size(pixelized , 1) , aux2(end) - m_x(end) )];

% making the size of the matrices consistent (Y direction)
aux1            = coef_g(1 , :);
aux2            = min(aux1(1) , m_y(1)) : max(aux1(end) , m_y(end));
coef_g_new      = [zeros(size(coef_g , 1)   , aux1(1) - aux2(1))   ,    coef_g   ,   zeros(size(coef_g , 1)   , aux2(end) - aux1(end))];
pixelized_new   = [zeros(aux2(end) - m_y(end) , size(pixelized_new , 2))   ;   pixelized_new  ;   zeros(m_y(1) - aux2(1) , size(pixelized_new , 2))];


% Defining "moments_dx_y"
moments_dx_y    = zeros(size(coef_dh , 1)-1  ,  size(coef_g , 1)-1);
for x_power = 0 : size(coef_dh , 1)-2
    
    aux         = sum( pixelized_new .* repmat(coef_dh_new(x_power+2 , :) , size(pixelized_new , 1) , 1)  ,  2 );
    
    for y_power = 0 : size(coef_g , 1)-2
        moments_dx_y(x_power+1 , y_power+1)     = coef_g_new(y_power+2 , end:-1:1) * aux;
    end
    
    
end








%%%%%%%
% Generating "moments_x_dy"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% making the size of the matrices consistent (X direction)
aux1            = coef_g(1 , :);
aux2            = min(aux1(1) , m_x(1)) : max(aux1(end) , m_x(end));
coef_g_new      = [zeros(size(coef_g , 1)    , aux1(1) - aux2(1))   ,    coef_g    ,   zeros(size(coef_g , 1)    , aux2(end) - aux1(end))];
pixelized_new   = [zeros(size(pixelized , 1) , m_x(1) - aux2(1) )   ,   pixelized  ,   zeros(size(pixelized , 1) , aux2(end) - m_x(end) )];

% making the size of the matrices consistent (Y direction)
aux1            = coef_dh(1 , :);
aux2            = min(aux1(1) , m_y(1)) : max(aux1(end) , m_y(end));
coef_dh_new     = [zeros(size(coef_dh , 1)   , aux1(1) - aux2(1))   ,    coef_dh   ,   zeros(size(coef_dh , 1)   , aux2(end) - aux1(end))];
pixelized_new   = [zeros(aux2(end) - m_y(end) , size(pixelized_new , 2))   ;   pixelized_new  ;   zeros(m_y(1) - aux2(1) , size(pixelized_new , 2))];


% Defining "moments_x_dy"
moments_x_dy    = zeros(size(coef_g , 1)-1  ,  size(coef_dh , 1)-1);
for x_power = 0 : size(coef_g , 1)-2
    
    aux         = sum( pixelized_new .* repmat(coef_g_new(x_power+2 , :) , size(pixelized_new , 1) , 1)  ,  2 );
    
    for y_power = 0 : size(coef_dh , 1)-2
        moments_x_dy(x_power+1 , y_power+1)     = coef_dh_new(y_power+2 , end:-1:1) * aux;
    end
    
    
end








%%%%%%%
% Generating "moments_x_y"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% making the size of the matrices consistent (X direction)
aux1            = coef_g(1 , :);
aux2            = min(aux1(1) , m_x(1)) : max(aux1(end) , m_x(end));
coef_g1_new     = [zeros(size(coef_g , 1)    , aux1(1) - aux2(1))   ,    coef_g    ,   zeros(size(coef_g , 1)    , aux2(end) - aux1(end))];
pixelized_new   = [zeros(size(pixelized , 1) , m_x(1) - aux2(1) )   ,   pixelized  ,   zeros(size(pixelized , 1) , aux2(end) - m_x(end) )];

% making the size of the matrices consistent (Y direction)
aux2            = min(aux1(1) , m_y(1)) : max(aux1(end) , m_y(end));
coef_g2_new     = [zeros(size(coef_g , 1)   , aux1(1) - aux2(1))   ,    coef_g   ,   zeros(size(coef_g , 1)   , aux2(end) - aux1(end))];
pixelized_new   = [zeros(aux2(end) - m_y(end) , size(pixelized_new , 2))   ;   pixelized_new  ;   zeros(m_y(1) - aux2(1) , size(pixelized_new , 2))];


% Defining "moments_x_y"
moments_x_y    = zeros(size(coef_g , 1)-1  ,  size(coef_g , 1)-1);
for x_power = 0 : size(coef_g , 1)-2
    
    aux         = sum( pixelized_new .* repmat(coef_g1_new(x_power+2 , :) , size(pixelized_new , 1) , 1)  ,  2 );
    
    for y_power = 0 : size(coef_g , 1)-2
        moments_x_y(x_power+1 , y_power+1)      = coef_g2_new(y_power+2 , end:-1:1) * aux;
    end
    
    
end