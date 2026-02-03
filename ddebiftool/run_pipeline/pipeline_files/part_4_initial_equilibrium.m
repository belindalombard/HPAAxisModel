% Follows part_3_functions.m in the pipeline
% File to define calculate the stability of the initial provided steady
% state, and attempt to correct it (get it closer to the actual steady
% state, if possible). 
% Output is two figures showing the initial steady state eigenvalues &
% corrected steady state eigenvalues

current_dir = pwd;
check_dir = 'pipeline_files';
if contains(current_dir, check_dir)
    cd('..')
end

%% Check that all the necessary parameters have been defined
required_params = {'funcs', 'stst'}; % list of parameters required for this script to run
check_parameters(required_params); % check that the required parameters are defined. 

%% Confirm that the output dir exists
if ~exist('results', 'dir')
   mkdir('results'); 
end
results_subdir = fullfile('results', filename);
if ~exist(results_subdir, 'dir')
   mkdir(results_subdir); 
end


%% Set some options that are necessary for the step
flag_newhheur=1; % Default choice

%% Define the method 
method=df_mthod(funcs,'stst',flag_newhheur);

method.continuation.step = 0.01;      
method.continuation.h_max = 0.01;    
method.newton.tol = 1e-6;             
method.newton.max_iter = 10;        

% Set the minimal real part for stability computations
method.stability.minimal_real_part=-1;

%% Correct the steady state
[stst, success] = p_correc(funcs, stst, [], [], method.point);

% If the steady state could not be found, indicate that. 
if ~success
    warning('Steady state correction failed. Proceeding with uncorrected state.');
end

%% Compute the stability
stst.stability=p_stabil(funcs,stst,method.stability);

% Display the calculated equilibrium (sanity check) 
fprintf('Equilibrium state vector: \n');
disp(stst.x); 

%% Plot the stability
figure(2); clf;
p_splot(stst); 
title('Stability of initial steady state');

saveas(gcf, fullfile(results_subdir, '1_initial_ss.png')); 

%% Ask for roots with more negative real part 
% Adjust minimal real part to identify more negative eigenvalues
method.stability.minimal_real_part=-3; 
stst.stability=p_stabil(funcs,stst,method.stability); % recompute stability:

%% Plot the stability of the new points
figure(3); clf;
p_splot(stst); 
title('Stability after correcting steady state'); % Title of plot

% Save plot
saveas(gcf, fullfile(results_subdir, '2_corrected_ss.png'));  % Save figure as PNG

% The next file to run is part_5_trivial_branch.m