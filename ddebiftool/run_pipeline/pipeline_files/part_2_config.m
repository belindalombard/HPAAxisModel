% Follows part_1_setup.m in the pipeline
% File to setup parameter values by reading values from a JSON
% configuration file. 

%% Read json file in (from config.json)
% Ensure that the path is correct relative to the current working directory.
% If you run into issues, run pwd to confirm the directory


current_dir = pwd;
check_dir = 'pipeline_files';
if contains(current_dir, check_dir)
    cd('..')
end

if ~exist('current_config_file', 'var') 
    disp('Running setup');
    run('pipeline_files/part_1_setup.m');
    current_config_file = 'config/config.json';
end

if ~exist('current_config_file_params', 'var') 
    current_config_file_params = 'config/config_parameters.json'
end

pwd

config = jsondecode(fileread(current_config_file));
config_params = jsondecode(fileread(current_config_file_params));

%% Extract necessary parameters from the config file 
init_parameters = config.initial_parameters; % contains all the initial params

% Confirm that parameters are in the correct order:
% Decide desired parameter order depending on whether one or two delays
if isfield(init_parameters, 'tau1') && isfield(init_parameters, 'tau2')
    desired_order = {'gamma_a','k_c','m_a','crh', ...
                     'gamma_c','alpha','m_c','k_a', ...
                     'tau1','tau2'};
else
    desired_order = {'gamma_a','k_c','m_a','crh', ...
                     'gamma_c','alpha','m_c','k_a','delay'};
end

parameters = struct();

for i = 1:numel(desired_order)
    parameters.(desired_order{i}) = init_parameters.(desired_order{i});
end

x_values = config.initial_values; % initial conditions for A and C
param_bif_1 = config.param_bif_1; % First bifurcation parameter (as a string)
param_bif_2 = config.param_bif_2; % Second bifurcation parameter (as a string)
filename = config.filename; 

index_bif_1 = find(strcmp(fieldnames(parameters), param_bif_1)); % Index of first bifurcation parameter
index_bif_2 = find(strcmp(fieldnames(parameters), param_bif_2)); % Index of second bifurcation parameter

if isfield(parameters, 'delay')
    % Single-delay case
    index_tau_c = find(strcmp(fieldnames(parameters), 'delay'));
elseif isfield(parameters, 'tau1') && isfield(parameters, 'tau2')
    % Two-delay case
    index_tau1 = find(strcmp(fieldnames(parameters), 'tau1'));
    index_tau2 = find(strcmp(fieldnames(parameters), 'tau2'));
end



% Parameters that define the bounds of the branches
n_steps = config.branch.n_steps;
min_bound_1 = config.branch.min_bound_1;
max_bound_1 = config.branch.max_bound_1; 
max_step_1 = config.branch.max_steps_1; 
max_iterations = config.branch.max_iterations; 
min_bound_2 = config.branch.min_bound_2;
max_bound_2 = config.branch.max_bound_2; 
max_step_2 = config.branch.max_steps_2; 

if ~exist('min_bound_1_plot', 'var')
    min_bound_1_plot = min_bound_1;
end

if ~isfield(config.branch, 'max_bound_1_plot')
    max_bound_1_plot = max_bound_1;
else
    max_bound_1_plot = config.branch.max_bound_1_plot;
end

if ~exist('min_bound_2_plot', 'var')
    min_bound_2_plot = min_bound_2;
end


if ~exist('max_bound_2_plot', 'var')
    max_bound_2_plot = max_bound_2;
end

if isfield(config_params, 'crh_parameters')
    crh_parameters = config_params.crh_parameters;
else
    crh_parameters = [];
end

if isfield(config_params, 'simulation_parameters')
    simulation_parameters = config_params.simulation_parameters;
else
    simulation_parameters = [];
end

%% Set JSON directory for reading parameter files
participant_json_dir = 'config/participants/'

%% Next file to run: part_3_functions.m