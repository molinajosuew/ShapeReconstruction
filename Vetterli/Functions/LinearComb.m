function [coefs , rel_error , max_error] = LinearComb(main_shape , n)
%
% [coefs , rel_error , max_error] = LinearComb(main_shape , n)
%
% This function finds the best linear combination of shifted splines of
% degree "n" that generate a given shape.
%
% "main_shape" is a 2*l matrix of the form
%               .-                                          -.
%               |     x_1          x_2      ...       x_l    |
%               |                                            |
%               | shape(x_1)   shape(x_2)   ...   shape(x_l) |
%               *-                                          -*
%              where the first row shall be an arithmetic progresion and 
%              the second row represents the values of the shape function 
%              at instances provided in the first row. 
%
%
% "n" is the degree of the basis splines used for representing the shape.
%
%
% "coefs" is a 2*L matrix that contains the coefficients in the linear
%         combination in its second row, while their indices is stored in
%         the first row:
%               .-                                       -.
%               |    m_1         m_2      ...      m_L    |
%               |                                         |
%               | coef(m_1)   coef(m_2)   ...   coef(m_L) |
%               *-                                       -*
%
% "rel_error" is the relative error of representing the "shape" using the
%             shifts of the spline; it is defined as the maximum absolute
%             error divided by the peak of the shape.
%
% "max_error" is the maximum error of representing the "shape" using the
%             shifts of the spline.


if (n < 0)||(n ~= round(n))
    error('!!! "n" should be a positive integer !!!')
end


fig_enable      = 1;        % enabling or disabling the plots
                            % 0: disable
                            % 1: enable


% primary parameters
dx              = main_shape(1,2) - main_shape(1,1);
m_range         = floor(main_shape(1,1) - (n+1)/2 - 2*n)  :  ceil(main_shape(1,end) + (n+1)/2 + 2*n);
x_range         = (m_range(1) - (n+1)/2) : dx : (m_range(end) + (n+1)/2) ;


% all the involved shifted splines
all_splines     = zeros(length(m_range) , length(x_range));
for m_ind = 1 : length(m_range)
    all_splines(m_ind , :)  = Centered_Spline(n , x_range - m_range(m_ind));
end


% pairwise inner-products between the shifted splines
spline_cros_cor = zeros(length(m_range));
for m_ind1 = 1 : length(m_range)
    for m_ind2 = 1 : length(m_range)
        spline_cros_cor(m_ind1 , m_ind2)    = (all_splines(m_ind1 , :) * all_splines(m_ind2 , :)') * dx;
    end
end


% extending the shape over 'x_range' by zero padding
shape_ext       = zeros(1 , length(x_range));
ind             = ( x_range >= main_shape(1,1)-1e-10 )&( x_range <= main_shape(1,end)+1e-10 );
shape_ext(ind)  = main_shape(2 , :);


% finding the inner-products between shifted splines and the shape
spline_shape_cor= zeros(length(m_range) , 1);
for m_ind = 1 : length(m_range)
    spline_shape_cor(m_ind)     = (all_splines(m_ind , :) * shape_ext') * dx;
end


% finding the optimal coefficients
opt_coef        = pinv(spline_cros_cor) * spline_shape_cor;


% forming the output
coefs           = [m_range ; opt_coef'];


% evaluating the error and plotting the approximated shape
if (nargout > 1)
    
    shape_app   = zeros(1 , length(x_range));
    for m_ind = 1 : length(m_range)
        shape_app   = shape_app + coefs(2 , m_ind) * all_splines(m_ind , :);
    end
    
    max_error   = max(abs(shape_ext-shape_app));
    rel_error   = max_error / max(shape_ext);
    
    if fig_enable == 1
        figure,
        subplot(2,1,1)
        plot(x_range , [shape_ext ; shape_app])
        legend('Original' , 'approximated')
        subplot(2,1,2)
        plot(x_range , shape_ext-shape_app)
        legend('error')
    end

    
end