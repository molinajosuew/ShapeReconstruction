function [pixelized , m_x , m_y] = pixelization(main_image , x_range , y_range , spline_degree)
%
% [pixelized , m_x , m_y] = pixelization(main_image , x_range , y_range , spline_degree)
%
% This function converts an image into pixels by applying a separable 2D
% PSF and sampling. The PSF is in fact the kronecker product of a B-spline 
% by itself.
%
% "main_image" is the original 2D image, possibly with high resolution.
%
% "x_range" and "y_range" are row vectors that reveal the grid over which
%           the original image is defined. Eventually, the sampling will 
%           take place at integer coordinates. Therefore, these two vectors
%           affect the number of pixels in the end.
%
% "spline_degree" is the degree of the spline used as the PSF of the
%                 imaging system responsible for generating the low 
%                 resolution image. It is assumed that the 2D PSF has a
%                 seprable form as a multiplication of a horizontal and a
%                 vertical B-spline of degree "spline_degree".
%
% "pixelized" is 2D matrix containing values between 0 and 1 that represent
%             the low resolution image. The number of elements depend on
%             "x_range" and "y_range".
%
% "m_x" and "m_y" store the X and Y location of the kernel centers for
%            sampling the algebraic cruve image.





% defining the sampling grid 
m_x         = ceil(min(x_range) - (spline_degree+1)/2) : floor(max(x_range) + (spline_degree+1)/2);
m_y         = ceil(min(y_range) - (spline_degree+1)/2) : floor(max(y_range) + (spline_degree+1)/2);


% shifting the x and y ranges according to the sampling grid
shifted_xx  = repmat(x_range , length(m_x) , 1) - repmat(m_x' , 1 , length(x_range));
shifted_yy  = repmat(y_range , length(m_y) , 1) - repmat(m_y' , 1 , length(y_range));

% separable implementation of the 2D kernel
aux_x       = Centered_Spline(spline_degree , shifted_xx);
aux_y       = Centered_Spline(spline_degree , shifted_yy);


% joint application of the kernel and the sampling
pixelized   = aux_y * main_image * aux_x' ...
              * (x_range(2) - x_range(1)) * (y_range(2) - y_range(1));



% plotting the low resolution image
if 0
    figure,
    imshow(imresize(pixelized , round(length(main_image) / length(pixelized)) , 'nearest') )
    
end