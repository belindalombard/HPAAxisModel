% Follows part_5_trivial_branch.m in the pipeline
% Locate the first Hopf point and continue the first Hopf bifurcation.
% Output is a bifurcation diagram of the hopf branch

current_dir = pwd;
check_dir = 'pipeline_files';
if contains(current_dir, check_dir)
    cd('..')
end
%% Check that all the necessary parameters have been defined
required_params = {'funcs', 'branch1', 'index_bif_1', 'index_bif_2', 'stst', 'config', 'n_steps'};
check_parameters(required_params); % Ensure all required parameters are present

%% Confirm that the output dir exists
if ~exist('results', 'dir')
   mkdir('results'); 
end
results_subdir = fullfile('results', filename);
if ~exist(results_subdir, 'dir')
   mkdir(results_subdir); 
end


%% Find the first hopf bifurcation point
% Find the last point in the branch from the previous step where we have
% Re(eigenvalue)>0. This is an index, of where this occurs. 
ind_hopf = find(arrayfun(@(x)real(x.stability.l0(1))>0, branch1.point), 1, 'last');

% Confirm that point has been found
if isempty(ind_hopf)
    error('No hopf point found in the branch. Check branch stability.');
end

% Convert the index to a hopf bifurcation point
hopf = p_tohopf(funcs, branch1.point(ind_hopf));

%% Correct and compute stability of the Hopf point
% Set up the method
method = df_mthod(funcs, 'hopf', flag_newhheur);
method.stability.minimal_real_part = -1;
method.continuation.step = 0.01;         
method.continuation.h_max = 0.01;        
method.newton.tol = 1e-6;                
method.newton.max_iter = 10;   

% "Correct" the hopf bifurcation point if needed
[hopf, success] = p_correc(funcs, hopf, index_bif_1, [], method.point); 

% Check if it succeeded
if ~success
    warning('Correction of hopf point failed.');
end

% Save first hopf bifurcation point
first_hopf = hopf; 

% Compute stability of hopf bifurcation point
hopf.stability = p_stabil(funcs, hopf, method.stability);

% Plot the stability of the hopf bifurcation point
figure(7); clf; 
p_splot(hopf); % plot stability of hopf bifurcation point
title('Stability of first hopf bifurcation'); 

% Save the figure
if ~exist('results', 'dir')
    mkdir('results');  % Create the 'results' directory if it does not exist
end
saveas(gcf, fullfile(results_subdir, '5_hopf_stability.pdf')); 

%% Initialize and continue the first Hopf bifurcation
% Extract bounds and steps from the config for branch 2
min_bound_1 = config.branch.min_bound_1;
min_bound_2 = config.branch.min_bound_2;
max_bound_1 = config.branch.max_bound_1;
max_bound_2 = config.branch.max_bound_2;
max_steps_1 = config.branch.max_steps_1;
max_steps_2 = config.branch.max_steps_2;

% Initialize the hopf branch using bif parameters 1 and 2
branch2 = df_brnch(funcs, [index_bif_1, index_bif_2], 'hopf'); 

% Set the bounds and step sizes for the bifurcation parameters
branch2.parameter.min_bound(1:2,:) = [[index_bif_1 min_bound_1]' [index_bif_2 min_bound_2]']';
branch2.parameter.max_bound(1:2,:) = [[index_bif_1 max_bound_1]' [index_bif_2 max_bound_2]']';
branch2.parameter.max_step(1:2,:) = [[index_bif_1 max_steps_1]' [index_bif_2 max_steps_2]']';

% Set the hopf point that we have found as the starting point of
% continuation
branch2.point = hopf;

% Slightly perturb the second bifurcation parameter for the next branch point
hopf.parameter(index_bif_2) = hopf.parameter(index_bif_2) + 0.01; 
[hopf, success] = p_correc(funcs, hopf, index_bif_1, [], method.point); 
if ~success
    warning('Correction hopf point failed.');
end

% Set the perturbed Hopf point as the second point of the branch
branch2.point(2) = hopf;

%% Plot the hopf bifurcation branch 
% Create the figure
figure(8); clf; set(gcf, 'color', 'w');

% Continue the hopf branch in one direction for n_steps
[branch2, s, f, r] = br_contn(funcs, branch2, n_steps); 

% Reverse the branch and continue in the opposite direction
branch2 = br_rvers(branch2); 
[branch2, s, f, r] = br_contn(funcs, branch2, n_steps); 

% Extract bifurcation parameter values from the branch points for plotting
for i = 1:length(branch2.point)
    x1(i) = branch2.point(i).parameter(index_bif_1);  % First bifurcation parameter
    y1(i) = branch2.point(i).parameter(index_bif_2);  % Second bifurcation parameter
end

% Plot the Hopf branch in blue
plot(x1, y1, '-b', 'LineWidth', 3); hold on;

% Add shading under the curve
fill([x1, max(x1)], [y1, max(y1)], [0 0.5 1], 'FaceAlpha', 0.3);

% Add dashed lines for x = 0 and y = 0 axes
plot([0, 0], [min_bound_2, max_bound_2], '--k', 'LineWidth', 1); 
plot([min_bound_1, max_bound_1], [0, 0], '--k', 'LineWidth', 1); 

% Add red line to show min CRH and max CRH
if exist('delay', 'var') && exist('crh_min', 'var') && exist('crh_max', 'var')
    plot([delay, delay], [crh_min, crh_max], '-r', 'LineWidth', 2); 
end

% Label the axes and set font sizes
xlabel(param_bif_1, 'FontSize', 24);
ylabel(param_bif_2, 'FontSize', 24);
title('Hopf bifurcation diagram');
axis([min_bound_1 max_bound_1 min_bound_2 max_bound_2]);
set(gca, 'fontsize', 18);


%% Save the figure
% Adjust properties 
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

% Configure save settings
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Save the figure 
savefile = '6_hopf_bifurcation_diagram';
saveas(fig, fullfile(results_subdir, [savefile, '.pdf']));


%% Next file to run: part_6b....m