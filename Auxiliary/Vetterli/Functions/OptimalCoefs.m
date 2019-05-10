function [b_i_m] = OptimalCoefs(Kernel , max_power , x_min , x_max)
%
% [b_i_m] = OptimalCoefs(Kernel , max_power , x_min , x_max)
% 
% This function evaluates the optimal coefficients in a linear combination
% of shifted kernels in order to generate polynomial modulations of a
% function. In simple words, we would like to find 'b(i,m)' not all zero 
% for which there exists a function 'g(x)' such that
%
%       d/(dx) (x^i g(x))  ~  \sum_{m} b(i,m) Kernel(x-m)
%
% where ~ shows almost equality.
%
%
% "Kernel" is a 2*l matrix that represents the kernel function. Its
%          strcuture is given by
%                        _                                           _ 
%                       |   x_1        x_2      x_3     ...    x_l    |
%              Kernel = |                                             |
%                       | Ker(x_1)  Ker(x_2)  Ker(x_3)  ...  Ker(x_l) |
%                       *-                                           -*
%          We assume that {x_i}s form an arithmetic progression. It is
%          further assumed that the kernel vanishes beyond the given range.
% 
%
% "max_power" is the maximum degree 'i' of the monomials multiplied by 
%             'g(x)'. Also, the minimum of 'i' is 1 and not 0. 
%
% "x_min" and "x_max" are the bounds on the interval for which we woul like
%         our almost equality to hold. Essentially, we do not care about 
%         the approximation quality outside of this range. These two values
%         should be integers.
%
%
% "b_i_m" is a 2D matrix containing the achieved optimal coefficients
%         b(i,m). This matrix includes "max_power" number of rows 
%         corresponding to i = 1 ... max_power. The number of columns is 
%         given by the number of integers between (not equal to)
%                   x_min - x_l     and     x_max - x_1
%                       


% primary parameters
dx              = Kernel(1,2) - Kernel(1,1);            % sampling resolution

x_range         = x_min : dx : x_max;                   % considered range of 'x' 

m_min           = floor(x_min - Kernel(1,end) + 1);     % minimum index 'm' for b(i,m)
m_max           = ceil(x_max - Kernel(1,1) - 1);        % maximum index 'm' for b(i,m)
m_range         = [m_min : m_max];                      % all values of 'm' for b(i,m)




% finding the inner product of the integrated kernels; i.e.,
%       \int_{x_min}^{x_max} x^i * int_ker_m1(x) * int_ker_m2(x) * dx
% where 
%           int_ker_m(x) = int_{0}^{x} Kernel(y-m) * dy
% and 'i' ranges from 0 to 2*(max_power - 1)
Ker_inner_prod  = zeros(length(m_range)  ,  length(m_range)  ,  2*max_power-1);

x_extended      = x_min + m_min : dx : x_max + m_max;
x_large         = x_min + 2*m_min : dx : x_max + 2*m_max;

Kernel_large    = zeros(size(x_large));
Kernel_large( ( x_large >= Kernel(1,1)-1e-10 )&( x_large <= Kernel(1,end)+1e-10 ) )     = Kernel(2,:);

Kernel_lar_int  = cumsum(Kernel_large) * dx;        % LSI integral of the kernel

if 1        % if we already know that the kernel is symmetric this part helps
    Kernel_lar_int  = 0.5*(Kernel_lar_int - Kernel_lar_int(end:-1:1) + 1);
end

for m_ind1 = 1 : length(m_range)
    aux         = (x_large + m_range(m_ind1) >= x_min+m_min-1e-10)&(x_large + m_range(m_ind1) <= x_max+m_max+1e-10);
    Ker_int_1   = Kernel_lar_int(aux) - Kernel_lar_int( abs(x_large + m_range(m_ind1)) < 1e-10);
    for m_ind2 = 1 : length(m_range)
        aux         = (x_large + m_range(m_ind2) >= x_min+m_min-1e-10)&(x_large + m_range(m_ind2) <= x_max+m_max+1e-10);
        Ker_int_2   = Kernel_lar_int(aux) - Kernel_lar_int( abs(x_large + m_range(m_ind2)) < 1e-10);
        
        x_power     = ones(size(x_extended));
        for i_ind = 1 : 2*max_power-1
            Ker_inner_prod(m_ind1 , m_ind2 , i_ind)     = sum(Ker_int_1 .* Ker_int_2 .* x_power) * dx;
            x_power         = x_power .* x_extended;
        end
        
    end
    
end




% Constructing the matrix of gradient equations
AA                  = zeros(max_power*length(m_range));
last_p              = 1 + (max_power-1) * length(m_range)  :  max_power * length(m_range);
for aa = 1 : max_power-1
    aux             = 1 + (aa-1) * length(m_range)  :  aa * length(m_range);
    AA(aux , aux)   = Ker_inner_prod(: , : , 2*(max_power-aa)+1);
    
    AA(aux , last_p)= - Ker_inner_prod(: , : , (max_power-aa)+1);
    
    AA(last_p , aux)= AA(aux , last_p).';
end
AA(last_p , last_p) = (max_power - 1) * Ker_inner_prod(: , : , 1);




% improving the condition number
AA_new              = zeros(size(AA));
for aa = 1 : length(AA)
    AA_new(aa , :)  = AA(aa , :) / norm(AA(aa , :));
end



% solving for the solution
[V , D]     = eig(AA_new);
b_i_m       = reshape(V(: , end) , length(m_range) , max_power)';