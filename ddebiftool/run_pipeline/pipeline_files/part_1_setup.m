% First file in the pipline
% A setup file to add the necessary paths, and clear the variables

%% Clear variables, format and add DDEBIFTOOL to the path.
clear all;             
format compact
% Add DDE-BIFTOOL to path — update this path to point to your local installation
if ispc
    addpath('C:\Users\belinda\Desktop\PhD\installation\DDE-Biftool\ddebiftool');
else
    % macOS / Linux: set DDEBIFTOOL_PATH env variable or edit the line below
    ddebiftool_path = getenv('DDEBIFTOOL_PATH');
    if isempty(ddebiftool_path)
        ddebiftool_path = fullfile(fileparts(fileparts(fileparts(pwd))), 'ddebiftool');
    end
    addpath(ddebiftool_path);
end

options = ddeset(...
    'MaxStep', 0.01, ...    
    'RelTol', 1e-6, ...     
    'AbsTol', 1e-9, ...    
    'InitialStep', 0.01);    

%% Solver defaults  
% (method is set later in part_3_functions_*.m once funcs is available)

%% Next file to run: part_2_config.m 

