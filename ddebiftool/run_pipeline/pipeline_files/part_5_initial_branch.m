% Follows part_4_initial_equilibrium.m in the pipeline
% Script to run intial branch of trivial equilibria. 

current_dir = pwd;
check_dir = 'pipeline_files';
if contains(current_dir, check_dir)
    cd('..')
end

%% Check that all the necessary parameters have been defined
required_params = {'n_steps', 'filename', 'min_bound_1', 'max_bound_1', 'max_step_1', 'max_iterations', 'index_bif_1', 'funcs', 'stst'}; % list of parameters required for this script to run
check_parameters(required_params); % check that the required parameters are defined. 

%% Confirm that the output dir exists
if ~exist('results', 'dir')
   mkdir('results'); 
end
results_subdir = fullfile('results', filename);
if ~exist(results_subdir, 'dir')
   mkdir(results_subdir); 
end


%% Initialise branch of trivial equilibria (the first branch) 
% Initialise the steady-state branch
branch1 = df_brnch(funcs,index_bif_1,'stst');

% Set bounds for the first bifurcation parameter
branch1.parameter.min_bound(1,:)=[index_bif_1 min_bound_1];
branch1.parameter.max_bound(1,:)=[index_bif_1 max_bound_1];
branch1.parameter.max_step(1,:)=[index_bif_1 max_step_1];

% use stst as a first branch point:
branch1.point=stst;

%%  Extend and continue branch of trivial equilibria
stst.parameter(index_bif_1)=stst.parameter(index_bif_1)+max_step_1;
[stst,success] = p_correc(funcs,stst,[],[],method.point);

if ~success
    warning('SS correction failed.');
end

% Use the corrected steady state as the second branch point
branch1.point(2)=stst;
branch1.method.continuation.plot=0; % Disable plotting during continuation

% Continue in one direction
% s is success, f is failed, and r is rejected
[branch1,s,f,r]=br_contn(funcs,branch1,max_iterations); 

%% Reverse the direction of the branch and continue
branch1=br_rvers(branch1); % Reverse the branch to go in opposite direction
[branch1,s,f,r]=br_contn(funcs,branch1,max_iterations);


%% Stability of branch of equilibria
% Set the minimal real part threshold for stability calculation.
branch1.method.stability.minimal_real_part=-2;

% Compute the eigenvals for each point along the branch
branch1=br_stabl(funcs,branch1,0,0);

% Obtain suitable scalar measures to plot stability along branch:
[xm,ym]=df_measr(1,branch1) % xm: bifurcation param values, ym: real part of eigenvalues

%% Plot stability along the branch of equilibria
figure(5); clf;
br_plot(branch1,xm,ym,'b'); 
ym.subfield='l0'; % Use real part of eigenvals
br_plot(branch1,xm,ym,'r'); % Line for eigenvalues (unstable)

%% Add additional things to plot to make it nicer
% Line at y=0, boundary for stability
plot([min_bound_1 max_bound_1], [0 0], '-.');
% Set x-axis limits. 
axis([min_bound_1 max_bound_1 -2 0.5]);
% Set x and y axis labels
xlabel(param_bif_1);  ylabel('\Re\lambda');
% Set plot title 
title('Branch stability.');

saveas(gcf, fullfile(results_subdir, '3_branch_stability.png'));  % Save the plot

%% Plot stability versus point number:
figure(6); clf;

% Plot the real part of eigenvalues versus point number along the branch
br_plot(branch1,[],ym,'b');
br_plot(branch1,[],ym,'b.');
% Set x and y label
xlabel('Point number along branch'); ylabel('\Re(\lambda)');
% Set plot ttile
title('Point number along the branch vs Eigenvalues');

% Save plot
saveas(gcf, fullfile(results_subdir, '4_stability_point_number.png'));  


% The next file to run is part_6_hopf_bifurcation
