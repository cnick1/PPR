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
%       Example 9: Chebychev discretization of Allen-Cahn equation
%
%
close all; clear; clc;
setKroneckerToolsPath
addpath('examples')

%% Example 11: Inverted Pendulum 
for nFterms = 3:9
        runExample11(nFterms, nFterms+1);
end

% Can run with acceleration using optional third argument
runExample11(9, 10, 1);

%% Example 9: Allen-Cahn equation, Dirichlet Boundary Conditions
% Allen-Cahn example where we wish to move the interface. This example
% features comparison with SDRE and the Bass 1966 approach which is like
% linear-polynomial regulation (sort-of). TT-HJB cannot be applied to this
% example because it is multi-input, and TT-HJB is only written for single
% input.
% close all

% Run cases with different diffusion coefficients
runExample9(45,4,10)
% runExample9(129,4,10)

% Run cases with different initial conditions
runExample9_differentICs(45, 4, 10)

%% Example 9: Allen-Cahn equation, Dirichlet Boundary Conditions
% clear;clc;close all;
runExample9_Neumann(14,6,0.2)

runExample9_Neumann(14,6,0.5)
