% Follows part_6_hopf_bifurcation.m in the pipeline
% Extra hopf bifurcation plots

current_dir = pwd;
check_dir = 'pipeline_files';
if contains(current_dir, check_dir)
    cd('..')
end

%% Check that all the necessary parameters have been defined
required_params = {'funcs', 'filename', 'branch1', 'index_bif_1', 'index_bif_2', 'stst', 'config', 'n_steps', 'x1'};
check_parameters(required_params); % Ensure all required parameters are present


%% Confirm that the output dir exists
if ~exist('results', 'dir')
   mkdir('results'); 
end

results_subdir = fullfile('results', filename);
if ~exist(results_subdir, 'dir')
   mkdir(results_subdir); 
end

%% Figure with y-axis on a log scale
figure(9); % Second figure for log scale
plot(x1, y1, '-b', 'LineWidth', 3); hold on;

% Add shading under the curve
fill([x1, max(x1)], [y1, max(y1)], [0 0.5 1], 'FaceAlpha', 0.3);

% Add dashed lines for x = 0 and y = 0 axes


plot([0, 0], [min_bound_2, max_bound_2], '--k', 'LineWidth', 1); 
plot([min_bound_1_plot, max_bound_1_plot], [0, 0], '--k', 'LineWidth', 1); 

% Set log scale for y-axis
set(gca, 'YScale', 'log');

y_ticks = [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4];
set(gca, 'YTick', y_ticks);

set(gca, 'YTickLabel', arrayfun(@(x) sprintf('10^{%d}', log10(x)), y_ticks, 'UniformOutput', false));


% Label the axes and set font sizes
xlabel(param_bif_1, 'FontSize', 24);
ylabel(param_bif_2, 'FontSize', 24);
title('Hopf bif (log scale)');


% axis([min_bound_1_plot max_bound_1_plot 0 max_bound_2]);

% Add horizontal line at y=0

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

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

savefile = '6b_hopf_bifurcation_diagram_log';
saveas(fig, fullfile(results_subdir, [savefile, '.png']));


% Next file to run is part_7_small_amplitude_orbit.m