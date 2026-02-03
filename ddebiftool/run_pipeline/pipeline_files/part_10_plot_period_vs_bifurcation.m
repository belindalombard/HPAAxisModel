% Follows part_9_induction_profiles.m in the pipeline
% Script to plot period dependency on both bifurcation parameters (index_bif_1 and index_bif_2).
% Output: Period vs bifurcation parameter plot for ACTH and CORT.

cd('..')

%% Check that all the necessary parameters have been defined
required_params = {'branch4', 'branch5', 'index_bif_1', 'index_bif_2'};
check_parameters(required_params); % Ensure all required parameters are loaded

%% Confirm that the output dir exists
if ~exist('results', 'dir')
   mkdir('results'); 
end
results_subdir = fullfile('results', filename);
if ~exist(results_subdir, 'dir')
   mkdir(results_subdir); 
end

%% Extract period data from branch4 and branch5 for plotting
% First bifurcation parameter (index_bif_1)
for i = 1:length(branch4.point)
    x3(i) = branch4.point(i).parameter(index_bif_1);  % First bifurcation parameter
    y3(i) = branch4.point(i).period;  % Period (T) for branch4
end

% Second bifurcation parameter (index_bif_2)
for i = 1:length(branch5.point)
    x4(i) = branch5.point(i).parameter(index_bif_2);  % Second bifurcation parameter
    y4(i) = branch5.point(i).period;  % Period (T) for branch5
end

%% Plot the period vs bifurcation parameters on the same axis
figure(13); clf; set(gcf, 'color', 'w'); % Create new figure with white background

% Plot period vs first bifurcation parameter (index_bif_1)
h1 = plot(x3, y3, 'LineWidth', 3, 'Color', [0, 0.4470, 0.7410]); hold on;

% Plot period vs second bifurcation parameter (index_bif_2)
h2 = plot(x4, y4, 'LineWidth', 3, 'Color', [0.8500, 0.3250, 0.0980]);

% Add legend and labels
legend([h1, h2], {param_bif_1, param_bif_2}, 'FontSize', 18);
xlabel('Bifurcation Parameter', 'FontSize', 24); % Label for x-axis
ylabel('Period (T)', 'FontSize', 24); % Label for y-axis
title('Period vs Bifurcation Parameters', 'FontSize', 24); % Add plot title

% Save the figure
saveas(gcf, fullfile(results_subdir, '9_period_vs_bifurcation.png')); % Save the figure
