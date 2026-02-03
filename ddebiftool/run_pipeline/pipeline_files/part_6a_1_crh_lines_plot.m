% Follows part_6_hopf_bifurcation.m in the pipeline
% Extra hopf bifurcation plots

current_dir = pwd;
check_dir = 'pipeline_files';
if contains(current_dir, check_dir)
    cd('..')
end
%% Check that all the necessary parameters have been defined
required_params = {'funcs', 'filename', 'config', 'figure_number'};
check_parameters(required_params); % Ensure all required parameters are present


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
xlabel('Delay ($\tau$)', 'Interpreter', 'latex');
ylabel('$log (\Phi)$ (constant)', 'Interpreter', 'latex');
title('', 'Interpreter', 'latex');


%% Add lines to figure
if exist('crh_data_per_participant', 'var') && ~isempty(crh_data_per_participant) && isstruct(crh_data_per_participant)
    % Sort by participant number
    participant_numbers = [crh_data_per_participant.participant_number];
    [~, sortIdx] = sort(participant_numbers);
    crh_data_per_participant = crh_data_per_participant(sortIdx);
    
    % 10 colours (for each participant) -- to avoid ugly blue/white that we
    % can't see. 
    custom_colors = [
        0.9, 0.4, 0.1;  
        0.3, 0.7, 0.3;  
        0.6, 0.3, 0.8; 
        0.8, 0.7, 0.2;  
        0.1, 0.6, 0.6; 
        0.7, 0.2, 0.4;  
        0.8, 0.5, 0.2; 
        0.6, 0.6, 0.2; 
        0.8, 0.4, 0.5;  
        0.6, 0.2, 0.2;  
    ];
    num_colors = size(custom_colors, 1);
    legend_entries = cell(1, length(crh_data_per_participant));
    plot_handles = gobjects(1, length(crh_data_per_participant)); 

    for i = 1:length(crh_data_per_participant)
        color_index = mod(i-1, num_colors) + 1;
        delay = crh_data_per_participant(i).delay;
        crh_min = crh_data_per_participant(i).min_crh;
        crh_max = crh_data_per_participant(i).max_crh;

        plot_handles(i) = plot([delay, delay], [crh_min, crh_max], 'Color', custom_colors(color_index, :), 'LineWidth', 2); 
        hold on;

        legend_entries{i} = ['Participant ', num2str(i)];
    end

    lgd = legend(plot_handles, legend_entries, 'Location', 'northeastoutside');
    hold off; 
else
    if exist('delay', 'var') && exist('crh_min', 'var') && exist('crh_max', 'var')
        plot([delay, delay], [crh_min, crh_max], '-r', 'LineWidth', 2); 
    end
end
xlim([10, 30]);
ylim([0 16000]);

savefile = '6_hopf_bifurcation_diagram_lines';


saveas(gcf, fullfile(results_subdir, [savefile, '.pdf']));
saveas(gcf, fullfile(results_subdir, [savefile, '.png']));

set(gca, 'XLabel', text('String', ''), 'YLabel', text('String', ''));
legend('off');
% Save the figure without labels
savefile_no_labels = '6_hopf_bifurcation_diagram_block_no_labels';
saveas(gcf, fullfile(results_subdir, [savefile_no_labels, '.pdf']));

