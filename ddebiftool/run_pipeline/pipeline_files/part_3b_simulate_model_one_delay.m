% Script to simulate DDE model

current_dir = pwd;
check_dir = 'pipeline_files';
if contains(current_dir, check_dir)
    cd('..')
end

if ~exist('results', 'dir')
   mkdir('results'); 
end

results_subdir = fullfile('results', filename);
if ~exist(results_subdir, 'dir')
   mkdir(results_subdir); 
end


addpath('functions')

n_days = 10;

%% Check that all the necessary parameters have been defined
required_params = {'parameters', 'index_tau_c', 'crh_parameters', 'simulation_parameters', 'filename'}; % list of parameters required for this script to run
check_parameters(required_params); % check that the required parameters are defined 

[fig, crh_min, crh_max, delay] = run_simulation(simulation_parameters, crh_parameters);

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Capture img as frame for future overlay
frame = getframe(gcf);    
img = frame.cdata;  

savefile = ['1_model_simulation' num2str(min_bound_1) '_' num2str(max_bound_1)];
saveas(fig, fullfile(results_subdir, [savefile, '.png']));

%% Next file to run: part_4_initial_equilibrium.m