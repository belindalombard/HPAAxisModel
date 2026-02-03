% Follows part_6_hopf_bifurcation.m in the pipeline
% Extra hopf bifurcation plots

current_dir = pwd;
check_dir = 'pipeline_files';
if contains(current_dir, check_dir)
    cd('..')
end
%% Check that all the necessary parameters have been defined
required_params = {'funcs', 'filename', 'config', 'img'};
check_parameters(required_params); % Ensure all required parameters are present


%% Confirm that the output dir exists
if ~exist('results', 'dir')
   mkdir('results'); 
end

results_subdir = fullfile('results', filename);
if ~exist(results_subdir, 'dir')
   mkdir(results_subdir); 
end

figure(figure_number);
axes('Position', [0.5, 0.5, 0.5, 0.5]);
imshow(img);

savefile = ['6_hopf_bifurcation_diagram_overlay'];
saveas(gcf, fullfile(results_subdir, [savefile, '.png']));
