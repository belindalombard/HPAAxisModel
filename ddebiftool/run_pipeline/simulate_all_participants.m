
% List all config files in the config directory
path_param_bifurcations = 'config/participants_bif';  % Get all JSON config files


for pn = 1:10
    % Set the current config file path as a global variable
    current_config_file = [path_param_bifurcations, '/Participant_bif_', int2str(pn), '.0.json' ];
    
    disp(['Running pipeline with config file: ', current_config_file]);

    if ~exist('current_config_file', 'var') == 1
        disp('Running setup');
        run('pipeline_files/part_1_setup.m')
        current_config_file = 'config/config.json'
    end

    current_config_file_params = ['config/participants/Participant_', num2str(pn), '.0.json'] 


    run('pipeline_files/part_2_config.m')

    run('pipeline_files/part_3_functions_one_delay.m')

    run('pipeline_files/part_3b_simulate_model_one_delay.m')
end 