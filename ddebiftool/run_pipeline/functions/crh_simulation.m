function crh = crh_simulation(t, crh_parameters)
    global crh_min crh_max;
    lambda_a = crh_parameters.lambda_a;
    sigma = crh_parameters.sigma;
    lambda_s = crh_parameters.lambda_s;
    T_c = crh_parameters.T_c;
    t_s = crh_parameters.t_s;
   
    crh = lambda_a * exp(lambda_s * cos(2 * pi * ((t - t_s) / T_c) + sigma * cos(2 * pi * ((t - t_s) / T_c))));
    
    crh_min = min(crh_min, crh);
    crh_max = max(crh_max, crh);
end