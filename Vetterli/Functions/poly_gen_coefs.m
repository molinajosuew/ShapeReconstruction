function [coef_g , error_g , coef_dg , error_dg] = poly_gen_coefs(spline_degree , g_degree , max_power , x_max)
%
% [coef_g , error_g , coef_dg , error_dg] = poly_gen_coefs(spline_degree , g_degree , max_power , x_max)
%
% The goal of this function is to represent 
%                   x^i * g(x)      and      d/dx (x^j * g(x))
% as a linear combination of shifted splines. Here, 'i' and 'j' range from 
% 0 to "max_power". Furthermore, 'g(x)' is itself a spline of degree
% "g_degree" with its supprt scaled to fit in the interval [-x_max , x_max].
% The splines used for representing 'x^i * g(x)' and 'd/dx (x^j * g(x))'
% are assumed to be of degree "spline_degree" in the standard form; i.e., 
% the spline shall have a support of length 'spline_degree+1'. The format 
% of the outputs is such that:
%
%          x^i * g(x) = \sum_{m}  coef_g(i+2 , m) * B(x - coef_g(1 , m))
%
%   d/dx (x^j * g(x)) = \sum_{m} coef_dg(j+2 , m) * B(x - coef_dg(1 , m))
%
% where 'B(x)' stands for the spline of order "spline_degree".
%
% "error_g" and "error_dg" are (max_power+1)*2 matrices whose first columns
% contain relative errors of representing 'x^i * g(x)' and 
% 'd/dx (x^j * g(x))' as a linear combination of shifted splines, 
% respectively. Similarly, the second columns contain maximum errors. 
%
% Note:
% "coef_g" and "coef_dg" have (max_power+2) rows while "error_g" and 
% "error_dg" have (max_power+1) rows. The reason is that the first rwo in
% both "coef_g" and "coef_dg" are reserved for the indices of the
% coefficients.





% defining the grid (by changing 'dx' one can modify the accuracy)
dx          = 1e-3;
x_range     = -x_max : dx : x_max;


% generating 'g(x)' and its derivative
scl         = (g_degree+1) / 2 / x_max;     % the 'x' scaling to fit the support
gx          = Centered_Spline(g_degree , x_range * scl);
dgx         = scl * (  Centered_Spline(g_degree-1 , x_range * scl + 0.5) ...
                     - Centered_Spline(g_degree-1 , x_range * scl - 0.5) );




% auxillary parameters
aux_dgx     = dgx;
aux_gx      = gx;


% defining the beginning of coef_dg
[coefs , rel_error , max_error] = LinearComb([x_range ; aux_dgx] , spline_degree);

coef_dg     = zeros(max_power+2 , size(coefs , 2));
coef_dg([1,2] , :)    = coefs;

error_dg    = zeros(max_power+1 , 2);
error_dg(1 , :)     = [rel_error , max_error];




% defining the beginning of coef_g
[coefs , rel_error , max_error] = LinearComb([x_range ; aux_gx] , spline_degree);

coef_g      = zeros(max_power+2 , size(coefs , 2));
coef_g([1,2] , :)       = coefs;

error_g     = zeros(max_power+1 , 2);
error_g(1 , :)      = [rel_error , max_error];




% sweepping over polynomial powers
for ii = 1 : max_power
    
    aux_dgx = aux_dgx .* x_range;
    [coefs , rel_error , max_error]     = LinearComb([x_range ; ii*aux_gx + aux_dgx] , spline_degree);
    
    coef_dg(ii+2 , :)   = coefs(2 , :);
    error_dg(ii+1 , :)  = [rel_error , max_error];
    
    
    
    aux_gx  = aux_gx  .* x_range;
    [coefs , rel_error , max_error]     = LinearComb([x_range ; aux_gx] , spline_degree);
    
    coef_g(ii+2 , :)    = coefs(2 , :);
    error_g(ii+1 , :)   = [rel_error , max_error];
end

