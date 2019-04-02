function cost = cost_dg(g0 , coef_dg , Kernel , lambda , max_power , x_min , x_max)
%
% cost = cost_dg(g0 , coef_dg , Kernel , lambda , max_power , x_min , x_max)
% 
% This function evaluates the cost incurred by representing the derivatives 
% of x^i g(x) by linear combinations of shifted kernels. 'g(.) is an
% unknown function. In ideal case, we would like to have 
%
%       d/(dx) (x^i g(x))  ~  \sum_{m} b(i,m) Kernel(x-m)
%
% where ~ shows almost equality and 'b(i,m)' are the coeeficients stored in
% "coef_dg". As we do not know 'g(.)', we only compare these equations with
% respect to each other:
%
%                                   cost = 
% sum_{i} ||\sum_{m} b(i,m) * K_m(x) - x * \sum_{m} b(i-1,m) * K_m(x)||_2^2
%
% where 'K_m(x)' stands for \int_{0}^{x} Kernel(y-m) * dy
%
%
%
%
%
% "g0" is ...
%
% "coef_dg" is a 2D matrix containing the coefficients b(i,m). This matrix 
%           includes "max_power+1" number of rows, the first of which
%           corresponds to values of 'm' and the following rows to monomial
%           degrees i = 1 ... max_power.
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
% "lambda" is a weight vector with length "max_power".
%
%
% "max_power" is the maximum degree 'i' of the monomials multiplied by 
%             'g(x)'. Also, the minimum of 'i' is 1 and not 0. 
%
% "x_min" and "x_max" are the bounds on the interval for which we woul like
%         our almost equality to hold. Essentially, we do not care about 
%         the approximation quality outside of this range. These two values
%         should be integers.                     


% primary parameters
dx              = Kernel(1,2) - Kernel(1,1);            % sampling resolution
dx              = round(dx * 1e7) / 1e7;

x_range         = x_min : dx : x_max;                   % considered range of 'x' 

m_min           = coef_dg(1,1);                         % minimum index 'm' for b(i,m)
m_max           = coef_dg(1,end);                       % maximum index 'm' for b(i,m)
m_range         = coef_dg(1,:);                         % all values of 'm' for b(i,m)
MMM             = length(m_range);




% evaluating the integrated kernels K_m(x) =\int_{0}^{x} Kernel(y-m) * dy
x_extended      = x_min + m_min : dx : x_max + m_max;
x_large         = x_min + 2*m_min : dx : x_max + 2*m_max;

Kernel_large    = zeros(size(x_large));
Kernel_large( ( x_large >= Kernel(1,1)-1e-10 )&( x_large <= Kernel(1,end)+1e-10 ) )     = Kernel(2,:);

Kernel_lar_int  = cumsum(Kernel_large) * dx;        % LSI integral of the kernel

Int_Ker_shift   = zeros(MMM  ,  length(x_extended));
for m_ind = 1 : MMM
    aux         = (x_large + m_range(m_ind) >= x_min+m_min-1e-10)&(x_large + m_range(m_ind) <= x_max+m_max+1e-10);
    Int_Ker_shift(m_ind , :)    = Kernel_lar_int(aux) - Kernel_lar_int( abs(x_large + m_range(m_ind)) < 1e-10);
end

aux             = (x_extended >= x_min-1e-10)&(x_extended <= x_max+1e-10);
Int_Ker_shift   = Int_Ker_shift(: , aux);





% implementing the cost function
cost            = lambda(1) * sum( ((g0 + coef_dg(2 , :) * Int_Ker_shift) .* x_range - coef_dg(3 , :) * Int_Ker_shift).^2 ) * dx;
for pp = 2 : max_power
    cost        = cost + ...
                  lambda(pp) * sum( ( (coef_dg(pp+1 , :) * Int_Ker_shift) .* x_range - coef_dg(pp+2 , :) * Int_Ker_shift).^2 ) * dx;
end

