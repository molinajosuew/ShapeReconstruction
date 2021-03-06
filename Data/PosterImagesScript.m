for N = 45 : 50
    %%%%%%%%%%%%%%
    % Parameters %
    %%%%%%%%%%%%%%
    
    n = 4;
    x_a = - 5;
    x_b = 5;
    y_a = - 5;
    y_b = 5;
    m_x = 2 ^ 7;
    m_y = 2 ^ 7;
    
    %%%%%%%%%%%%
    % Computed %
    %%%%%%%%%%%%
    
    L = x_b;
    x_n = m_x * L + 1;
    y_n = m_y * L + 1;
    
    %%%%%%%%%%
    % Images %
    %%%%%%%%%%  
    
    I = GetImageOfPower(GetRandomPower(n, x_a, x_b, y_a, y_b, false), x_a, x_b, y_a, y_b, x_n, y_n);
%     I = GetImageOfNonSeparable(GetRandomNonSeparable(n, L), L, x_n, y_n);
    
    %%%%%%%%%%
    % Saving %
    %%%%%%%%%%
    
    path = "/home/wjmolina/Documents/GitHub/ShapeReconstruction/Data/PosterImages/";
    mkdir(sprintf(path + "%i", N));
    imwrite(ind2rgb(im2uint8(mat2gray(I)), parula(256)), path + N + "/originalImage.png");
    
    for psnr = [- 1, 30, 50, 100]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Down-Sample and Reconstructions %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [R_V, ~, D] = VetterliReconstruction(I, x_a, x_b, y_a, y_b, x_n, y_n, psnr);
        
        best = - 1;
        for j = 1 : 100
%             R_B_candidate = GetImageOfNonSeparable(NonSeparableReconstruction(I, n, L, psnr), L, x_n, y_n);
            [CP, ~] = PowerReconstruction(I, n, x_a, x_b, y_a, y_b, psnr);
            R_B_candidate = GetImageOfPower(CP, x_a, x_b, y_a, y_b, x_n, y_n);
            if SorensenDiceCoefficient(I, R_B_candidate) > best
                best = SorensenDiceCoefficient(I, R_B_candidate);
                R_B = R_B_candidate;
            end
        end
        
        %%%%%%%%%%
        % Saving %
        %%%%%%%%%%
        
        if psnr == - 1
            psnr = "Noiseless";
        end
        
        imwrite(ind2rgb(im2uint8(mat2gray(D)), parula(256)), path + N + "/sampledImage.png");
        imwrite(ind2rgb(im2uint8(mat2gray(abs(I - R_B))), parula(256)), path + N + "/bernsteinReconstruction" + psnr + "PSNR.png");
        imwrite(ind2rgb(im2uint8(mat2gray(abs(I - R_V))), parula(256)), path + N + "/vetterliReconstruction" + psnr + "PSNR.png");
    end
end