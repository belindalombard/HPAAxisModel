cd('..')

%% Construction and continuation of branch5 (periodic solutions with respect to the second bifurcation parameter)
% Initialize an empty branch for periodic solutions along the second bifurcation parameter
branch5 = df_brnch(funcs, index_bif_2, 'psol'); % Initialize branch5 for bifurcation parameter 2

% Set bounds for continuation along the second bifurcation parameter
branch5.parameter.min_bound(1,:) = [index_bif_2 min_bound_2];
branch5.parameter.max_bound(1,:) = [index_bif_2 max_bound_2];
branch5.parameter.max_step(1,:) = [index_bif_2 max_steps_2];

% Use degenerate periodic solution and corrected periodic solution as starting points
deg_psol = p_topsol(funcs, first_hopf, 0, degree, intervals); % Degenerate periodic solution
deg_psol.mesh = []; % Clear mesh
branch5.point = deg_psol; % Set degenerate solution as the first point

psol.mesh = []; % Clear mesh for the periodic solution
branch5.point(2) = psol; % Set corrected periodic solution as the second point

%% Continue branch5 for periodic solutions along the second bifurcation parameter
figure(11); clf; % Create new figure for continuation of branch5
[branch5, s, f, r] = br_contn(funcs, branch5, 70); % Continue branch5 for periodic solutions


