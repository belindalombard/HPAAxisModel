function [fig, crh_min, crh_max, delay] = run_simulation(simulation_parameters, crh_parameters)
    n_days = 10;
    delay = simulation_parameters.delay;
    tspan = [1 1440 * n_days];
    y0 = @(t) [1; 0];

    % Initialize CRH tracking variables
    global crh_min crh_max;
    crh_min = Inf;
    crh_max = -Inf;
    crh_values = [];

    % Generate CRH values for the simulation period
    for t = 1:1:1441
        crh_value = crh_simulation(t, crh_parameters);
        crh_values(end + 1) = crh_value;
    end
    

    sol = dde23(@(t, y, ylag) dde_system(t, y, ylag, simulation_parameters, crh_parameters), delay, y0, tspan);
        

    
    t_last_day = (n_days * 1440 - 1440):1:(n_days * 1440);
    y_last_day = deval(sol, t_last_day);

    fig = figure(1);

    % Plot A on the left y-axis
    yyaxis left
    plot(t_last_day, y_last_day(1, :), 'b-', 'LineWidth', 1.5); 
    ylabel('ACTH');
    ylim([min(y_last_day(1, :)) - 1, max(y_last_day(1, :)) + 1]);

    yyaxis right
    plot(t_last_day, y_last_day(2, :), 'r-', 'LineWidth', 1.5); 
    hold on;
    plot(t_last_day, crh_values, 'g--', 'LineWidth', 1.5);
    ylabel('CORT');
    ylim([min([y_last_day(2, :) crh_values]) - 1, max([y_last_day(2, :) crh_values]) + 1]); 
    hold off;

    % Legend, labels, and title
    legend('ACTH', 'CORT', 'CRH', 'Location', 'best');
    xlabel('Time (minutes)');
    title('Simulation of A and C');
    grid on;
end