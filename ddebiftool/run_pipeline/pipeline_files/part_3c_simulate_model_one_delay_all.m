% Script to simulate DDE model for all participants.
% Script saves the minimum and maximum CRH values for each participant, so
% it can be used. 

%% Check that all the necessary parameters have been defined
required_params = {'parameters', 'index_tau_c', 'filename', 'participant_json_dir'}; % list of parameters required for this script to run
check_parameters(required_params); % check that the required parameters are defined 

%% Confirm the current directory, and change if necessary
current_dir = pwd;
check_dir = 'pipeline_files';
if contains(current_dir, check_dir)
    cd('..')
end

%% Define results directory and create if not exists
if ~exist('results', 'dir')
   mkdir('results'); 
end

results_subdir = fullfile('results', filename);
if ~exist(results_subdir, 'dir')
   mkdir(results_subdir); 
end

%% Add path for functions (Crh and dde function) 
addpath('functions')

%% Get JSON files for each paarticipant
json_files = dir(fullfile(participant_json_dir, '*.json'));

%% Set empty structure to save CRH data for each participant
crh_data_per_participant = struct('participant_number', {}, 'min_crh', {}, 'max_crh', {}, 'delay', {});

%% Set how long simulation should run
n_days = 10;

for k = 1:length(json_files)
    % Load JSON file
    json_file = fullfile(participant_json_dir, json_files(k).name);
    data = jsondecode(fileread(json_file));

    results_subdir = fullfile('results', filename);
    if ~exist(results_subdir, 'dir')
        mkdir(results_subdir); 
    end


    % Keep track of min and max CRH values
    global crh_min crh_max;
    crh_min = Inf;
    crh_max = -Inf;
    crh_values = [];

    participant_number = data.participant_number; 

   
    % Create figure
    figure(1);

    [fig, crh_min, crh_max, delay] = run_simulation(data.simulation_parameters, data.crh_parameters);
    

    % Save min and max CRH
    crh_data_per_participant(k).participant_number = participant_number;
    crh_data_per_participant(k).min_crh = crh_min;
    crh_data_per_participant(k).max_crh = crh_max;
    crh_data_per_participant(k).delay = data.simulation_parameters.delay;

    title(['Participant ' num2str(participant_number) ]);

    savefile = ['1_model_simulation' num2str(participant_number)];
    saveas(fig, fullfile(results_subdir, [savefile, '.png']));


end

%% Next file to run: part_4_initial_equilibrium.m