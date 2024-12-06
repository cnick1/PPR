%  A script to run the examples in SIAM SISC paper
%  The examples are the following:
%   
%       Example 10: 1D ODE toy example inspired by Jeff's 2018 paper. The code
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
addpath('utils')

%% Example 10: 1D Toy Example
% Prints out the closed-loop dynamics of LQR and PPR controllers
runExample10();

%% Example 11: Inverted Pendulum 
% Produces the value function plots, residual plots, and closed-loop phase
% portrait (by calling external python script)
for nFterms = 1:2:9
    runExample11(nFterms, nFterms+1);
end
% Can run with acceleration using optional third argument to do reduction
runExample11(9, 10, 1);

% Can play around with RofX and QofX for pendulum example, didn't have much success
% for nFterms = 1:2:9
%     runExample11_RofX(nFterms, nFterms+1);
% end

%% Example 9: Allen-Cahn equation, Dirichlet Boundary Conditions
% Allen-Cahn example where we wish to move the interface. This example
% features comparison with SDRE and the Bass 1966 approach which is like
% linear-polynomial regulation (sort-of). TT-HJB cannot be applied to this
% example (TT-HJB is only written for single input).
close all; clear; clc

% Run cases with different diffusion coefficients
runExample9(45,4,10)
% runExample9(129,4,10)

% Run cases with different initial conditions
runExample9_differentICs(45, 4, 10)

%% Example 9: Allen-Cahn equation, Neumann Boundary Conditions
% clear;clc;close all;
runExample9_Neumann(14,6,0.5)

%% Example 22: 4D Pendulum Cart Example
runExample22();