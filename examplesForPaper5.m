%  A script to run the examples in PPR paper for CDC2024
%  The examples are the following:
%
%       Example 7: 3D aircraft stall model for testing stabilization. From
%       Garrard 1977, also studied in Almubarak 2019 ACC paper and others.
%
%       Example 9: Chebychev discretization of Allen-Cahn equation
% 
close all; clear; clc;
setKroneckerToolsPath
addpath('examples')

%% Example 7: 3D aircraft stall model
runExample7()

%% Example 9: Allen-Cahn equation
runExample9(129,2)
close all
runExample9(129,4)
close all
