function [z_total,w_vector,n_begin,n_end,R_lb,R_ub,tD_lb,tD_ub,lD_lb,lD_ub]...
    = DDT_LoadSiNWData(filename_data, c_Li)
    

    % This function returns the selected part (at a specific concentration)
    % of the EIS data set.
    % Also returns appropriate hyperparameters for the quadratic
    % programming and the cross-validation

    data = load(filename_data);    
       
    
    if c_Li == 954 %1
        % Selected data
        w_vector= data.w_0954; % frequency data (vector)
        z_total = data.z_0954; % impedance data (vector, complex)
        n_begin = 38; % transition data point, the first data point to include in DDT analysis
        n_end = 54;  % the last data point to include in DDT analysis
        % hyperparameters
        R_lb = 4.34; % real-axis intercept of diffusion impedance, lower bound in fitting
        R_ub = 4.4; % real-axis intercept of diffusion impedance, uppder bound in fitting  
        tD_lb = 0.5; % time scale to include in DDT analysis, lower limit, log scale   
        tD_ub = 4; % time scale to include in DDT analysis, upper limit, log scale  
        lD_lb = -2; % regularization parameter (lambda), lower limit, log scale
        lD_ub = -1; % regularization parameter (lambda), upper limit, log scale  
         
    elseif c_Li==1274 %2
        % Selected data
        w_vector= data.w_1274;
        z_total = data.z_1274;
        n_begin = 39; % 37 -> 39
        n_end = 57; 
        % hyperparameters
        R_lb = 4.45;  
        R_ub = 4.47; 
        tD_lb = 0.5;  
        tD_ub = 4.5;    
        lD_lb = -2;  
        lD_ub = -0.5;   
    
    elseif c_Li== 2385 %3
        % Selected data
        w_vector= data.w_2385;
        z_total = data.z_2385;
        n_begin = 40; 
        n_end = 54;  
        % hyperparameters
        R_lb = 4.62;  
        R_ub = 4.73; 
        tD_lb = -1; 
        tD_ub = 6;  
        lD_lb = -2 ; % lD = -1 gives you a single peak
        lD_ub = 0; 
        
    elseif c_Li==2705 %4
        % Selected data
        w_vector= data.w_2705;
        z_total = data.z_2705;
        n_begin = 41; 
        n_end = 57;  
        % hyperparameters
        R_lb = 5.1;  
        R_ub = 5.3; 
        tD_lb = -2; 
        tD_ub = 5; 
        lD_lb = -3 ;
        lD_ub = 0; 
    
    end
    

end