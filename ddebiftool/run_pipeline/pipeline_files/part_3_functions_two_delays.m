% Follows part_2_config.m in the pipeline
% File to define the function, as well as the necessary parameters of the
% functions (initial parameters, delays)
% configuration file.
% Output is a stst object with the defined x and params :)

current_dir = pwd;
check_dir = 'pipeline_files';
if contains(current_dir, check_dir)
    cd('..')
end

addpath('functions')

%% Check that all the necessary parameters have been defined
% Now we need two delay indices (tau1: C->A feedback; tau2: A->C adrenal)
required_params = {'parameters', 'index_tau1', 'index_tau2', 'x_values'}; 
check_parameters(required_params); % ensure all required parameters are present

%% Define the system functions
% xx is ACTH and CORT.
%  - First dimension: state variable (1 = ACTH, 2 = CORT)
%  - Second dimension:
%       xx(:,1)   -> current state at time t
%       xx(:,2)   -> state at time t - tau1  (used for C feedback onto A)
%       xx(:,3)   -> state at time t - tau2  (used for A driving C in adrenal)
%
% Parameter vector 'par' (by index) is assumed as:
%   par(1) = gamma_a      (ACTH decay)
%   par(2) = K_c          (C half-saturation in feedback on A)
%   par(3) = m_a          (Hill exponent for C feedback on A)
%   par(4) = Phi          (constant CRH drive used for continuation)
%   par(5) = gamma_c      (CORT decay)
%   par(6) = alpha        (max adrenal production gain)
%   par(7) = m_c          (Hill exponent for A->C activation)
%   par(8) = K_a          (A half-saturation in adrenal activation)
%   par(index_tau1) = tau1 (C->A feedback delay)
%   par(index_tau2) = tau2 (A->C adrenal delay)
%
% NOTE: For DDE-BIFTOOL we keep Phi constant (no explicit t-dependence).

AC_sys_RHS = @(xx, par) [ ...
    -par(1) * xx(1,1) + (par(2)^(par(3)) * par(4)) ./ (par(2)^(par(3)) + (xx(2,2)).^(par(3))); ...
    -par(5) * xx(2,1) + par(6) * ( (xx(1,3)).^(par(7)) ) ./ ( par(8)^(par(7)) + (xx(1,3)).^(par(7)) ) ...
];

%% Define the delays
% Return the parameter indices for tau1 and tau2
AC_tau = @() [index_tau1, index_tau2];

%% Create a map of the defined parameters (handy utility, unchanged)
param_map = containers.Map(fieldnames(parameters), struct2cell(parameters));

%% Define structure of functions
funcs = set_funcs( ...
    'sys_rhs', AC_sys_RHS, ...
    'sys_tau', AC_tau);

%% Create a SS object
% This is the "initial guess" for a steady state solution
stst.kind      = 'stst';
stst.parameter = struct2array(parameters);     % row vector
stst.x         = reshape(struct2array(x_values), [], 1);  % column vector

%% Next file to run: part_4_initial_equilibrium.m
