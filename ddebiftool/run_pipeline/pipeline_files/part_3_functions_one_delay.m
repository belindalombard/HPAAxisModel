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
required_params = {'parameters', 'index_tau_c'}; % list of parameters required for this script to run
check_parameters(required_params); % check that the required parameters are defined 

%% Define the system functions
% This is the system. xx is ACTH and CORT. In the first dimension, it
% represents the current state of ACTH and CORT respectively. In the 2nd
% dimension, it represents the state of ACTH and CORT delayed with the
% first delay parameter (defined in the next step), etc. E.g. xx(2,3) is
% CORT tau_c steps back. 
AC_sys_RHS = @(xx, par) [
    -par(1) * xx(1, 1) + (par(2)^(par(3)) * par(4)) / (par(2)^(par(3)) + xx(2, 2)^(par(3)));...  
    -par(5) * xx(2, 1) + par(6) * (xx(1, 1)^par(7)) / (par(8)^(par(7)) + (xx(1, 1)^par(7))) 
];


%% Define the delays
% Here, 9 and 10 represents the 9th and 10th positions of the parameters
% object, and consequently tau_a and tau_c. 
% TODO - INDEXES DYNAMIC 
AC_tau=@()[index_tau_c];

%% Create a map of the defined parameters
% This map object maps parameter names to parameter values, allowing easy
% access to the parameters later on.
param_map = containers.Map(fieldnames(parameters), struct2cell(parameters));

%% Define structure of functions
funcs = set_funcs(...
    'sys_rhs', AC_sys_RHS,...
    'sys_tau', AC_tau);  

%% Create a SS object
% This is the "initial guess" for a steady state solution
stst.kind='stst';
stst.parameter = struct2array(parameters); % It has to be a row vector (',') 
stst.x = reshape(struct2array(x_values), [], 1); % It has to be a column vector (';')

%% Next file to run: part_4_initial_equilibrium.m