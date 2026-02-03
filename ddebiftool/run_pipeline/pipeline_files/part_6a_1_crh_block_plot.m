%% Confirm that the output dir exists
if ~exist('results', 'dir')
   mkdir('results'); 
end

results_subdir = fullfile('results', filename);
if ~exist(results_subdir, 'dir')
   mkdir(results_subdir); 
end

fig = figure(figure_number);


%% Set scale to normal if not set
if ~exist('scale', 'var')
    scale = 'normal';
end

%% Determine the overall range for delay and CRH
min_delay = inf;
max_delay = -inf;
min_crh = inf;
max_crh = -inf;

if exist('crh_data_per_participant', 'var') && ~isempty(crh_data_per_participant) && isstruct(crh_data_per_participant)
    % Iterate over each participant to find the min/max values for delay and CRH
    for i = 1:length(crh_data_per_participant)
        delay = crh_data_per_participant(i).delay;
        crh_min = crh_data_per_participant(i).min_crh;
        crh_max = crh_data_per_participant(i).max_crh;
        
        % Update overall min/max for delay and CRH
        min_delay = min(min_delay, delay);
        max_delay = max(max_delay, delay);
        min_crh = min(min_crh, crh_min);
        max_crh = max(max_crh, crh_max);
    end
end


% Plot the block covering the entire range
%hold on;
%h = fill([min_delay, max_delay, max_delay, min_delay], ...
%     [min_crh, min_crh, max_crh, max_crh], 'r', 'FaceAlpha', 0.5);
%hold off;

% Set axis labels using LaTeX interpreter for the first figure
xlabel('Delay ($\tau$) (minutes)', 'Interpreter', 'latex');
ylabel('CRH $(\Phi)$', 'Interpreter', 'latex');
title('', 'Interpreter', 'latex');

% xmax = max(max_delay, 10);   
xlim([15, 30]);



hold on;
x_dots = 22.14 * ones(1,3);
y_dots = [0.1, 10, 10000];
p = plot(x_dots, y_dots, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);  % black filled circles
hold off;

ylim([0 16000]);


%legend([h p], {'Individual trajectories region', 'Constant CRH simulations'}, 'Location', 'northoutside');

% Save the figure with labels
savefile = '6_hopf_bifurcation_diagram_block';
print(fullfile(results_subdir, savefile), '-dpdf', '-r300');  % Save the figure with labels

% Remove the x and y axis labels
% set(gca, 'XLabel', text('String', ''), 'YLabel', text('String', ''));
%legend('off');
%% Save the figure without labels
%savefile_no_labels = '6_hopf_bifurcation_diagram_block_no_labels';
% print(fullfile(results_subdir, savefile_no_labels), '-dpdf', '-r300');  