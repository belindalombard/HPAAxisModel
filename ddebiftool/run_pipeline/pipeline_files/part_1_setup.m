% First file in the pipline
% A setup file to add the necessary paths, and clear the variables

%% Clear variables, format and add DDEBIFTOOL to the path.
clear all;             
format compact
addpath('C:\Users\belinda\Desktop\PhD\installation\DDE-Biftool\ddebiftool');   

options = ddeset(...
    'MaxStep', 0.01, ...    
    'RelTol', 1e-6, ...     
    'AbsTol', 1e-9, ...    
    'InitialStep', 0.01);    

method = df_mthod('stst');

method.stability.minimal_real_part = -1e-5;   
method.point.extra_condition = 1e-6;          
method.continuation.step = 0.01;              
method.continuation.h_max = 0.01;             
method.newton.tol = 1e-6;                    
method.newton.max_iter = 10;                  

%% Next file to run: part_2_config.m 

