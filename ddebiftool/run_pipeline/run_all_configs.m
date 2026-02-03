% List all config files in the 'config' directory
config_files = dir('config/matlab/*.json');  % Get all JSON config files

for i = 1:length(config_files)
    % Set the current config file path as a global variable
    current_config_file = fullfile(config_files(i).folder, config_files(i).name);
    
    disp(['Running pipeline with config file: ', current_config_file]);
    
    run_steps; 
end
