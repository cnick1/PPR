%  A script to run the examples in SIAM SISC paper
%  The examples are the following:
%   
%       Example 10: 1D ODE toy example inspired by Jeff's 2018 paper. The code
%       here is just to compute energy functions and plot controllers
%       working/not working. In the write-up I hardcoded the energy
%       function values, so there is nothing to export with this example.
%
%       Example 26: Nondimensionalized 2D inverted pendulum. 
%
%       Example 9: Chebychev discretization of Allen-Cahn equation
%
%
close all; clear; clc;
setKroneckerToolsPath
addpath('examples')
addpath('utils')

[~,sys] = memory; RAM = sys.PhysicalMemory.Total/1e9; 

%% Example 10: 1D Toy Example with Transcritical Bifurcation
% Prints out the closed-loop dynamics of LQR and PPR controllers
runExample10();

%% Example 26: 2D Inverted Pendulum 
% Produces the value function plots, residual plots, and closed-loop phase
% portraits (by calling external python script)
for nFterms = 1:2:9
    runExample26(nFterms+1);
end

% Prints out the table of control costs
runExample26_controlCosts();

%% Example 9: Allen-Cahn equation, Neumann Boundary Conditions
% Case where we wish to stabilize the origin and compare with other
% methods, most notably TT-HJB. 
runExample9_Neumann(14,6,0.5)

%% Example 9: Allen-Cahn equation, Dirichlet Boundary Conditions
% Case where we wish to move the interface. TT-HJB cannot be applied to 
% this example because i) it only is written for single input systems, and
% ii) it is too high dimensional

if RAM < 20
    n = 45;          % for running locally     with  16 GB RAM
else 
    n = 129;         % for running on Poseidon with 512 GB RAM
end

% Run cases with different diffusion coefficients
runExample9(n,4,10)            

% Run cases with different initial conditions
runExample9_differentICs(n,4,10)

%% Example 22: 4D Pendulum Cart Example
% runExample22();