% Follows part_7_small_amplitude_orbits.m in the pipeline
% Construct and continue a branch of periodic orbits.
% Output: Periodic orbit branch, amplitude plots for ACTH and CORT

cd('..')

%% Check that all the necessary parameters have been defined
required_params = {'funcs', 'index_bif_1', 'filename', 'min_bound_1', 'max_bound_1', 'max_steps_1', 'psol', 'first_hopf'};
check_parameters(required_params); % Ensure necessary parameters are in scope

%% Confirm that the output dir exists
if ~exist('results', 'dir')
   mkdir('results'); 
end
results_subdir = fullfile('results', filename);
if ~exist(results_subdir, 'dir')
   mkdir(results_subdir); 
end


%% Construction and continuation of the periodic orbit branch
% Set up the branch
branch4 = df_brnch(funcs, index_bif_1, 'psol'); 

% Set the bounds for continuation along the bifurcation parameter
branch4.parameter.min_bound(1,:) = [index_bif_1 min_bound_1];
branch4.parameter.max_bound(1,:) = [index_bif_1 max_bound_1];
branch4.parameter.max_step(1,:) = [index_bif_1 max_steps_1];

% Create a periodic solution with zero amplitude at the hopf point
deg_psol = p_topsol(funcs, first_hopf, 0, degree, intervals); 
deg_psol.mesh = []; % for continuation

% Set solutions as the first two points of the branch
branch4.point = deg_psol; 
psol.mesh = []; 
branch4.point(2) = psol; % Corrected solution as second point


%% Continue the branch
% Create figure for branch continuation
figure(10); clf; set(gcf, 'color', 'w'); 
[branch4, s, f, r] = br_contn(funcs, branch4, n_steps); % Compute periodic solutions 

% Plot amplitude of periodic solutions for ACTH and CORT
for i = 1:length(branch4.point)
    x2(i) = branch4.point(i).parameter(index_bif_1); 
    y21(i) = max(branch4.point(i).profile(1,:)) - min(branch4.point(i).profile(1,:));  % Amplitude of ACTH
    y22(i) = max(branch4.point(i).profile(2,:)) - min(branch4.point(i).profile(2,:));  % Amplitude of CORT
end

h1 = plot(x2, y21, 'LineWidth', 3, 'Color', [0, 0.4470, 0.7410]); hold on;
h2 = plot(x2, y22, 'LineWidth', 3, 'Color', [0.8500, 0.3250, 0.0980]); hold on;

% Add figureinfo
legend([h1 h2], {'ACTH', 'CORT'}, 'FontSize', 18, 'Position', [0.15 0.7 0.2 0.2]);
legend('boxoff'); % Remove legend box
xlabel(param_bif_1, 'FontSize', 24); % Label for bifurcation parameter
ylabel('Amplitude', 'FontSize', 24); % Label for amplitude
title('Amplitude of Periodic Orbits for ACTH and CORT');

% Save amplitude plot
saveas(gcf, fullfile(results_subdir, '7_periodic_orbit_amplitude.png'));

% Next file to run: part_9_period_vs_bifurcation.m