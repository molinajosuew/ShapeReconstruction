function [coef_dg , error_dg] = OptimalCoefs_dg(Kernel , lambda , max_power , x_min , x_max)
%
% [coef_dg , error_dg] = OptimalCoefs_dg(Kernel , lambda , max_power , x_min , x_max)
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
%
%
% "coef_dg" is a 2D matrix containing the achieved optimal coefficients
%           b(i,m). This matrix includes "max_power" number of rows 
%           corresponding to i = 1 ... max_power. The number of columns is 
%           given by the number of integers between (not equal to)
%                     x_min - x_l     and     x_max - x_1
%                       
% "error_dg" is a (max_power+1) * 2 matrix that contains the max and
%            relative absolue errors in representing x^i g(x) using the
%            shifted kernels. The first column contains the relative
%            errors. Also, the i'th row corresponds to  x^(i-1) g(x).


% primary parameters
dx              = Kernel(1,2) - Kernel(1,1);            % sampling resolution
dx              = round(dx * 1e7) / 1e7;

x_range         = x_min : dx : x_max;                   % considered range of 'x' 

m_min           = floor(x_min - Kernel(1,end) + 1);     % minimum index 'm' for b(i,m)
m_max           = ceil(x_max - Kernel(1,1) - 1);        % maximum index 'm' for b(i,m)
m_range         = [m_min : m_max];                      % all values of 'm' for b(i,m)
MMM             = length(m_range);




% extending the grid
x_extended      = x_min + m_min : dx : x_max + m_max;
x_large         = x_min + 2*m_min : dx : x_max + 2*m_max;

Kernel_large    = zeros(size(x_large));
Kernel_large( ( x_large >= Kernel(1,1)-1e-10 )&( x_large <= Kernel(1,end)+1e-10 ) )     = Kernel(2,:);



% LSI integral (from -Inf) of the extended kernel
Kernel_lar_int  = cumsum(Kernel_large) * dx;        % LSI integral of the kernel

if 1        % if we already know that the kernel is symmetric this part helps
    Kernel_lar_int  = 0.5*(Kernel_lar_int - Kernel_lar_int(end:-1:1) + 1);
end



% shifted integrated kernels
Int_Ker_shift   = zeros(MMM  ,  length(x_extended));
for m_ind = 1 : MMM
    aux         = (x_large + m_range(m_ind) >= x_min+m_min-1e-10)&(x_large + m_range(m_ind) <= x_max+m_max+1e-10);
    Int_Ker_shift(m_ind , :)    = Kernel_lar_int(aux) - Kernel_lar_int( abs(x_large + m_range(m_ind)) < 1e-10);
end






% finding the inner product of the integrated kernels; i.e.,
%       \int_{x_min}^{x_max} x^i * int_ker_m1(x) * int_ker_m2(x) * dx
% where 
%           int_ker_m(x) = int_{0}^{x} Kernel(y-m) * dy
% and 'i' ranges from 0 to 2
Ker_Ker_prod    = zeros(MMM  ,  MMM  ,  3);

indicator       = (x_extended >= x_min-1e-10)&(x_extended <= x_max+1e-10);
%indicator       = ones(size(x_extended));

for m_ind1 = 1 : MMM
    for m_ind2 = 1 : MMM
        
        x_power     = ones(size(x_extended));
        for i_ind = 1 : 3
            Ker_Ker_prod(m_ind1 , m_ind2 , i_ind)       = sum(   Int_Ker_shift(m_ind1 , :) .* Int_Ker_shift(m_ind2 , :) ...
                                                              .* x_power .* indicator) * dx;
            x_power         = x_power .* x_extended;
        end
        
    end
    
end





% finding the inner product of the integrated kernels with monomials; i.e.,
%       \int_{x_min}^{x_max} x^i * int_ker_m(x) * dx
% where 
%           int_ker_m(x) = int_{0}^{x} Kernel(y-m) * dy
% and 'i' ranges from 1 to 2
Ker_x_prod      = zeros(2 , MMM);

for m_ind = 1 : MMM
    
    x_power     = ones(size(x_extended));
    for i_ind = 1 : 2
        x_power = x_power .* x_extended;
        Ker_x_prod(i_ind , m_ind)       = sum(Int_Ker_shift(m_ind , :) .* x_power .* indicator) * dx;
    end
        
end









%%%%%%%%%%%%%%%%%
% Constructing the Hessian matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AA                  = zeros((max_power+1) * MMM + 1);


% the exceptional first element
AA(1 , 1)           = sum( (x_extended .* indicator).^2 ) * dx;


% the first row and column
aux_ind             = 2 : MMM+1;

AA(1 , aux_ind)             =  lambda(1) * Ker_x_prod(2 , :);
AA(1 , aux_ind + MMM)       = -lambda(1) * Ker_x_prod(1 , :);

AA(aux_ind , 1)             = AA(1 , aux_ind)';
AA(aux_ind + MMM , 1)       = AA(1 , aux_ind + MMM)';


% the first MMM rows after the 1st row
AA(aux_ind , aux_ind)       = lambda(1) * Ker_Ker_prod(: , : , 3);
AA(aux_ind , aux_ind + MMM) = -lambda(1) * Ker_Ker_prod(: , : , 2);


% middle part of the matrix
for aa = 2 : max_power
    aux_ind         = 2 + MMM*(aa-1) : 1 + MMM*aa;
    
    AA(aux_ind , aux_ind - MMM) = -lambda(aa-1) * Ker_Ker_prod(: , : , 2);
    AA(aux_ind , aux_ind)       =  lambda(aa-1) * Ker_Ker_prod(: , : , 1) + lambda(aa) * Ker_Ker_prod(: , : , 3);
    AA(aux_ind , aux_ind + MMM) = -lambda(aa)   * Ker_Ker_prod(: , : , 2);
    
end


% last MMM rows
aux_ind                     = 2 + MMM*max_power : 1 + MMM*(max_power+1);
AA(aux_ind , aux_ind - MMM) = -lambda(max_power) * Ker_Ker_prod(: , : , 2);
AA(aux_ind , aux_ind)       =  lambda(max_power)   * Ker_Ker_prod(: , : , 1);

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













%%%%%%%%%%%%%%%%%
% minimizing the constrained QP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
norm(AA - AA')
min(eig(AA))



% denoising 'AA' by removing non-zero eigen-values
[V , D]     = eig(AA);
BB          = AA - min(min(D)) * eye(size(AA));



% parameters for quadratic programming
H           = BB;
f           = zeros(length(BB) , 1);

Aeq         = [1 zeros(1 , length(BB)-1)];
beq         = [1];

Aineq       = -[zeros(1 , length(x_extended)) ; ...
                Int_Ker_shift ; ...
                zeros(max_power * MMM , length(x_extended))]';
bineq       = ones(length(x_extended) , 1);

x0          = V(: , 1);

options     = optimset('MaxIter' , 1e6);


% QP for obtaining the coefficients
raw_coefs   = quadprog(H , f , Aineq , bineq , Aeq , beq , [] , [] , x0 , options);


% reshaping the coeffieicnts to form the final output
coef_dg     = [m_range ; reshape(raw_coefs(2:end) , MMM , max_power+1)'];

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%
% Checking the error 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_dg    = zeros(max_power+1 , 2);

% reconstructing the original 'g(.)' function
ggg         = coef_dg(2 , :) * Int_Ker_shift + raw_coefs(1);

% 'x' indices from the extended grid that coincide with the given range
aux         = (x_extended >= x_min-1e-10)&(x_extended <= x_max+1e-10);


% plotting 'g(.)'
figure(1)
plot(x_range , ggg(aux))

% plotting 'x^i * g(x)' and finding the error based on 'g(.)'
for pp = 1 : max_power
    ggg     = x_extended .* ggg;
    series  = coef_dg(pp+2 , :) * Int_Ker_shift;
    
    f1      = ggg(aux);
    f2      = series(aux);

    figure(pp+1)
    plot(x_range , [f1 ; f2])

    error_dg(pp+1 , :)      = [max(abs(f1-f2))/max(abs(f1))  max(abs(f1-f2))];
end