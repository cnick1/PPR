%  A script to run the examples in CONTROLS paper
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
%
%       Pendulum/acrobat, etc. ?
%       Example 5: a unicycle quadratic-bilinear system.
%
%       Example 6: a nonlinear (von Karman) Euler-Bernoulli beam with cable
%       actuation that provides state dependent inputs. The drift and input
%       dynamics are up-to cubic. The model can be made arbitrarily large.
%
%       Example 8: a nonlinear heat equation. The dynamics are cubic. The
%       model can be made arbitrarily large. According to Embree 2019,
%       etc., the open-loop linearization should be STABLE, but due to
%       non-normal A matrix the nonlinear system is actually UNSTABLE. 
%
close all; clear; clc;
setKroneckerToolsPath
addpath('examples')

exportData = true;

%% Example 10: 1D ODE
% runExample10_regionOfAccuracy_Lyapunov()

%% Example 11: 2D inverted pendulum
% runExample11(exportPlotData, nFterms, degree)
% runExample11_plotEnergyFunctions(exportPlotData, nFterms, degree, eta)

% for degree = 2:2:10
    % nFterms = degree - 1;
        % Phase portraits
%         runExample11(exportData, nFterms, degree)
%         runExample11_computeCost(false, nFterms, degree)
        % Value Functions
        % runExample11_plotEnergyFunctions(exportData, nFterms, degree)
% end


%% Example 7: 3D aircraft stall model
% runExample7(exportData)

%% Example 9: Allen-Cahn equation
runExample9()

runExample9(45,2)
runExample9(45,4)

runExample9(33,2)
runExample9(33,4)
runExample9(33,6)

%%
runExample9(129,2)
runExample9(129,4)

%% Other not good examples
%% Example 5: 3D unicycle, not locally stabilizable for small eps
% for degree = 2:2:10
%         runExample5(false, degree-1, degree, 0.001, 1)
% end


