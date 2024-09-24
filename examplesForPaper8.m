%  A script to run the examples in CONTROLS+MOR paper
%  The examples are the following:
%   
%       Example 10: 1D ODE toy example from Jeff's 2018 paper. The code
%       here is just to compute energy functions and plot controllers
%       working/not working. In the write-up I hardcoded the energy
%       function values, so there is nothing to export with this example.
%
%       Example 11: 2D inverted pendulum. We can choose different degree
%       approximations to the dynamics, different degree energy functions
%       to compute, etc.
%
%       Example 7: 3D aircraft stall model for testing stabilization. From
%       Garrard 1977, also studied in Almubarak 2019 ACC paper and others.
% 
%       Example 9: Chebychev discretization of Allen-Cahn equation
%
%
close all; clear; clc;
setKroneckerToolsPath
addpath('examples')
addpath('utils')
addpath('../NLbalancing')

exportData = true;

%% Example 9: Allen-Cahn equation
% close all
% runExample9_mor(45,45,2)
% runExample9_mor(45,10,2)
% 
runExample9_mor(45,45,4)
runExample9_mor(45,10,4)

runExample9_mor_u(129,10,2)
runExample9_mor_u(129,10,4)
runExample9_mor_u(129,10,6)

%% Example 11: Inverted Pendulum 
runExample11_mor(false, 7, 6,2)
runExample11_mor(false, 7, 6,1)

runExample11_mor_plotEnergyFunctions(exportData, 2, 2,2)
runExample11_mor_plotEnergyFunctions(exportData, 11, 10,2)
runExample11_mor_plotEnergyFunctions(exportData, 11, 10,1)

%% Example 11: Inverted Pendulum, Radial Basis Functions
for epsilon = 1:0.25:10
    runExample11_plotEnergyFunctions_RBF(epsilon)
end

%%
epsilons = [0.0001 0.01 0.05 0.1 0.25 0.5 1 2 5 10 100];
for epsilon = epsilons
    runExample11_plotEnergyFunctions_RBF_sigmoids(epsilon)
end