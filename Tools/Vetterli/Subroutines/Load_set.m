% This file is to read the right set of moment generating coefficients

if coef_setting_method == 1     % the coefficients are set based on 
                                % "spline_degree" and "sum_degree" automatically 
    
    maxPower    = ceil(sum_degree * 1.5);
    
    switch spline_degree
        case 2
            load(['coefs_n_2_m_6_maxPower_' , num2str(maxPower) , '_xMax_14.mat']);
            
        case 4
            load(['coefs_n_4_m_7_maxPower_' , num2str(maxPower) , '_xMax_16.mat']);
            
        case 6
            load(['coefs_n_6_m_10_maxPower_' , num2str(maxPower) , '_xMax_16.5.mat']);
            
        case 8
            load(['coefs_n_8_m_10_maxPower_' , num2str(maxPower) , '_xMax_16.5.mat']);
    end
    
    
    
else        % user-defined coefficients
    
    
    filename    = 'coefs_n_8_m_8_maxPower_8_xMax_13.5_shft_n_2.mat';
    load(filename);
    
    
end