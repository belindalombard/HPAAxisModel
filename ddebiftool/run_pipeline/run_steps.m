% Navigate to the correct directory
current_dir = pwd;
desired_dir = 'run_pipeline';
check_pipeline_dir = 'attempt_1';
if ~contains(current_dir, desired_dir)
    if contains(current_dir, check_pipeline_dir)
        cd(['../', desired_dir]);
    else 
        cd(desired_dir);
    end
else
    disp('Already in the desired directory.');
end


if ~exist('current_config_file', 'var') == 1
    disp('Running setup');
    run('pipeline/part_1_setup.m')
    current_config_file = 'config/config.json'
end


run('pipeline_files/part_2_config.m')

run('pipeline_files/part_3_functions_one_delay')

% Look at stability of equilibrium 
run('pipeline_files/part_4_initial_equilibrium.m')

% Look at initial branch
run('pipeline_files/part_5_initial_branch.m')


% Perform hopf bifurcation
run('pipeline_files/part_6_hopf_bifurcation.m')

% Plot hopf bifurcation on log scale
run('pipeline_files/part_6b_hopf_bifurcation_logscale_plot.m')

run('pipeline_files/part_7_small_amplitude_orbit')

run('pipeline_files/part_8_periodic_orbit')

run('pipeline_files/part_8b_second_continuation')

run('pipeline_files/part_9_period_vs_bifurcation')