% Navigate to the correct directory
%current_dir = pwd;
%desired_dir = 'run_pipeline';
%check_pipeline_dir = 'attempt_1';
%if ~contains(current_dir, desired_dir)
%    if contains(current_dir, check_pipeline_dir)
%        cd(['../', desired_dir]);
%    else 
%        cd(desired_dir);
%    end
%else
%    disp('Already in the desired directory.');
%end


if ~exist('current_config_file', 'var') == 1
    disp('Running setup');
    run('pipeline_files/part_1_setup.m')
    current_config_file = 'config/config.json'
end

if ~exist('current_config_file_params', 'var') == 1
    current_config_file_params = 'config/config_parameters.json'
end

run('pipeline_files/part_2_config.m')

run('pipeline_files/part_3_functions_one_delay.m')

if ~exist('crh_data_per_participant', 'var') || isempty(crh_data_per_participant) || ~isstruct(crh_data_per_participant)
    run('pipeline_files/part_3c_simulate_model_one_delay_all.m')
end 

% Look at stability of equilibrium 
run('pipeline_files/part_4_initial_equilibrium.m')

% Look at initial branch
run('pipeline_files/part_5_initial_branch.m')


% Perform hopf bifurcation
run('pipeline_files/part_6_hopf_bifurcation.m')

% Plot hopf bifurcation on log scale
run('pipeline_files/part_6b_hopf_bifurcation_logscale_plot.m')

figure_number = 8; 
run('pipeline_files/part_6a_1_crh_block_plot.m');

figure_number = 9;
run('pipeline_files/part_6a_1_crh_block_plot.m'); 
% run('pipeline_files/part_6c_overlay_plot.m')
