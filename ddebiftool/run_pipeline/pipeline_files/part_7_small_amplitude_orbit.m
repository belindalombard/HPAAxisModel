% Follows part_6_hopf_bifurcation.m in the pipeline
% Construct a small-amplitude periodic orbit near a hopf bifurcation. The
% output is a periodic solution correction

cd('..')

%% Check that all the necessary parameters have been defined
required_params = {'funcs', 'title', 'first_hopf', 'index_bif_1', 'min_bound_1', 'max_bound_1', 'max_steps_1'};
check_parameters(required_params); % Ensure necessary parameters are loaded


%% Constructing an initial small-amplitude orbit near a Hopf bifurcation
% The periodic solution is initialized with small amplitude (1e-2)
intervals = 18;
degree = 3;
[psol, stepcond] = p_topsol(funcs, first_hopf, 1e-2, degree, intervals); % Create periodic solution

% Correct the periodic solution guess
% Define method for periodic solutions
method = df_mthod(funcs, 'psol'); 
% Correct solution
[psol, success] = p_correc(funcs, psol, 10, stepcond, method.point);
if ~success
    warning('Correction of periodic solution failed.');
end

% Next file to run: part_8_periodic_orbit_branch.m
