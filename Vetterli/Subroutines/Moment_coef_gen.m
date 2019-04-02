% This file generates spline coefficients for generating 'x^i * g(x)' and 
% its derivative. The spline is assumed to be of degree "n" and centered 
% around the origin. Function 'g(.)' is a scaled version of another 
% centered spline (degree "m"). The scaling is such that the support spans
% from "-x_max" up to "x_max". It is also possible to include a time shift
% be setting the value in "shft". The range of 'i' considered for 
% 'x^i * g(x)' is from 0 to "max_power".


clc
clear all
close all







current_path    = pwd;
root_path       = current_path(1 : end-11);
addpath([root_path , 'Functions'])







% main parameters
n               = 8;
m               = 10;
x_max           = 16.5;
shft            = 0;

max_power       = 8;



% dependent parameters
dx              = 1e-3;
x_range         = -x_max : dx : x_max;
        
        
        




% generating 'g(x)' and its derivative
scl             = (m+1) / 2 / x_max;     % the 'x' scaling to fit the support
gx              = Centered_Spline(m , x_range * scl);
dgx             = scl * (  Centered_Spline(m-1 , x_range * scl + 0.5) ...
                         - Centered_Spline(m-1 , x_range * scl - 0.5) );

% generating the coefficients
[coef_g , error_g , coef_dg , error_dg] = spline_gen_coefs(n , [x_range+shft ; gx ; dgx] , max_power);





% saving the coefficients
if shft > 0
    shft_text   = ['_shft_p_' , num2str(shft)];
elseif shft<0
    shft_text   = ['_shft_n_' , num2str(abs(shft))];
else
    shft_text   = [];
end
save([root_path , 'Moment_Coefs/coefs_n_' , num2str(n) , '_m_' , num2str(m) , '_maxPower_' , num2str(max_power) , '_xMax_' , num2str(x_max) , shft_text , '.mat'] , 'coef_g' , 'error_g' , 'coef_dg' , 'error_dg')
