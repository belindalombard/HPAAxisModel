% Follows part_8_periodic_orbit_branch.m in the pipeline
% Script to plot the induction profiles of ACTH and CORT at the final point of the branch.
% Output: Induction profiles for ACTH and CORT.

cd('..')


%% Check that all the necessary parameters have been defined
required_params = {'branch4', 'index_bif_1'};
check_parameters(required_params); % Ensure all required parameters are loaded

%% Confirm that the output dir exists
if ~exist('results', 'dir')
   mkdir('results'); 
end
results_subdir = fullfile('results', filename);
if ~exist(results_subdir, 'dir')
   mkdir(results_subdir); 
end

%% Plot the profiles (induction) of ACTH and CORT at the end of branch4
figure(12); clf;  set(gcf, 'color', 'w'); % Create a new figure with white background

% Extract the final point from branch4 for plotting
bif_param_1_end = length(branch4.point); 

% Define the time grid for plotting the profiles
t = linspace(0, 1, length(branch4.point(bif_param_1_end).profile));

% Plot the profile of ACTH (first state variable)
g1 = plot(t, branch4.point(bif_param_1_end).profile(1,:), 'LineWidth', 3, 'Color', [0, 0.4470, 0.7410]); hold on;

% Plot the profile of CORT (second state variable)
g2 = plot(t, branch4.point(bif_param_1_end).profile(2,:), 'LineWidth', 3, 'Color', [0.8500, 0.3250, 0.0980]); hold on;

% Add legend for ACTH and CORT profiles
legend([g1 g2], {'ACTH', 'CORT'}, 'FontSize', 18, 'Position', [0.75 0.7 0.2 0.2]);
legend('boxoff'); % Remove the box around the legend

% Add axis labels
xlabel('Time / T', 'FontSize', 24); % Time (normalized by period T)
ylabel('Induction', 'FontSize', 24); % Induction level for ACTH and CORT
title('Induction Profiles for ACTH and CORT'); % Add a title

%% Automatically scale the y-axis to fit both profiles
ylim([min(min(branch4.point(bif_param_1_end).profile)) max(max(branch4.point(bif_param_1_end).profile))]);

% Set x-axis limits to range from 0 to 1 for normalized time
xlim([0 1]);

% Adjust figure layout to remove extra white space
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

%% Save the induction profile plot
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(gcf, fullfile(results_subdir, '8_induction_profiles_ACTH_CORT_fixed.png')); % Save the figure
