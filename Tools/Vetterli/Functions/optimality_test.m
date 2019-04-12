clc
clear all
close all





% main parameters
Spl_deg         = 8;

x_max           = 16.5;
x_min           = -16.5;

max_power       = 8;

lambda          = ones(1,max_power);



% dependent parameters
dx              = 1e-2;
x_range         = x_min : dx : x_max;
        
Kernel          = [x_range ; Centered_Spline(Spl_deg , x_range)];



% cost handle
c_dg            = @(coefs) cost_dg(1 , coefs , Kernel , lambda , max_power , x_min , x_max);
c_g             = @(coefs) cost_g(coefs , Kernel , lambda , max_power , x_min , x_max);





% loading a data set
%load('coefs_n_6_m_6_maxPower_6_xMax_10.5.mat')
load('coefs_n_8_m_10_maxPower_8_xMax_16.5.mat')

coef_gL     = coef_g(2:end , :) / coef_g(2 , coef_g(1 , :) == 0);
coef_gL     = [coef_g(1 , :)  ;  coef_gL];

coef_dgL    = coef_dg(2:end , :) / Centered_Spline(10 , 0);
coef_dgL    = [coef_dg(1 , :) ; coef_dgL];





%%

% optimization for 'g'
[coef_gO , error_gO]    = OptimalCoefs_g(Kernel , lambda , max_power , x_min , x_max);
coef_gO_modif   = coef_gO(2:end , :) / coef_gO(2 , coef_gO(1 , :) == 0);
coef_gO_modif   = [coef_gO(1,:) ; coef_gO_modif];

% cost for this choice
c_g(coef_gO_modif)

% cost for the loaded data
c_g(coef_gL)




%%

% optimization for 'dg'
[coef_dgO , error_dgO]  = OptimalCoefs_dg(Kernel , lambda , max_power , x_min , x_max);

% cost for this choice
c_dg(coef_dgO)

% cost for the loaded data
c_dg(coef_dgL)
























