%  A script to run the examples in SIAM SISC paper
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

exportData = true;


%% Example 11: Inverted Pendulum 
for nFterms = 3:9
        runExample11(nFterms, nFterms+1);
end

% Can run with acceleration using optional third argument
runExample11(9, 10, 1);

%% Example 9: Allen-Cahn equation
% close all
runExample9(45,4)
runExample9(45,4,10)
runExample9(129,4,10)

for x0 = -.75:.25:.75
    runExample9_differentICs(33, 4, x0)
end


