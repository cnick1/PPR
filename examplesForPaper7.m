%  A script to run the examples in PPR Pade/Chisholm/Canterbury paper
%  The examples are the following:
%
%       Example 11: 2D inverted pendulum. We can choose different degree
%       approximations to the dynamics, different degree energy functions
%       to compute, etc.
%
close all; clear; clc;
setKroneckerToolsPath
addpath('examples')

exportData = false;

%% Example 11: 2D inverted pendulum
% Ell = 3
runExample11_chisholmApproximant(3,2,1) 
runExample11_chisholmApproximant(3,3,2)

% Ell = 5
runExample11_chisholmApproximant(5,4,1)
runExample11_chisholmApproximant(5,5,2)

runExample11_chisholmApproximant(5,2,1)
runExample11_chisholmApproximant(5,3,2)
runExample11_chisholmApproximant(5,4,3)
runExample11_chisholmApproximant(5,5,4) %

% Ell = 7
runExample11_chisholmApproximant(7,6,1)
runExample11_chisholmApproximant(7,7,2)

runExample11_chisholmApproximant(7,4,1)
runExample11_chisholmApproximant(7,5,2)
runExample11_chisholmApproximant(7,6,3) %
runExample11_chisholmApproximant(7,7,4)

runExample11_chisholmApproximant(7,2,1)
runExample11_chisholmApproximant(7,3,2)
runExample11_chisholmApproximant(7,4,3)
runExample11_chisholmApproximant(7,5,4)
runExample11_chisholmApproximant(7,6,5)
runExample11_chisholmApproximant(7,7,6)




