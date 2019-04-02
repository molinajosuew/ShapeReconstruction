function beta_n_x = Centered_Spline(n , x)
%
% beta_n_x = Centered_Spline(n , x)
%
% This function generates the values of the "n"th degree spline at points
% included in "x". Thus, "beta_n_x" has the same size as "x".


if (n < 0)||(n ~= round(n))
    error('!!! "n" should be a positive integer !!!')
end

% employing the symmetry of the splines
x               = -abs(x);

aux1            = x + (n+1)/2;
aux1(aux1 < 0)  = 0;
aux2            = 1;

beta_n_x        = aux1 .^ n;
for kk = 1 : n+1
    aux1            = x + (n+1)/2 - kk;
    aux1(aux1 < 0)  = 0;
    aux2            = aux2 * (n+2-kk) / kk;
    
    beta_n_x    = beta_n_x + (-1)^kk * aux2 * (aux1 .^ n);
end

%beta_n_x(x >= (n+1)/2)  = 0;
%beta_n_x(beta_n_x < 0)  = 0;
beta_n_x                = beta_n_x / factorial(n);