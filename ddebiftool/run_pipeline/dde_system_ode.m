function dydt = dde_system_ode(t, y, history_func, delay, parameters, crh_parameters)
    A = y(1);
    C = y(2);

    % Get delayed value of C
    C_delay1 = history_func(t - delay, 2); 

    gamma_a = parameters.gamma_a;
    K_c = parameters.k_c;
    m_a = parameters.m_a;
    gamma_c = parameters.gamma_c;
    alpha = parameters.alpha;
    m_c = parameters.m_c;
    K_a = parameters.k_a;

    % Define system dynamics
    dAdt = -gamma_a * A + ((K_c^m_a) * crh_simulation(t, crh_parameters)) / (K_c^m_a + C_delay1^m_a);
    dCdt = -gamma_c * C + alpha * ((A^m_c) / (K_a^m_c + A^m_c));  
    
    dydt = [dAdt; dCdt];
end
