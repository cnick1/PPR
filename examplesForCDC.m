%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Main script to run the examples for CDC %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Author: Nicholas Corbin, UCSD           %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Usage: Run the entire script to generate the examples for the CDC 2024
%         Milan paper [1]. 
% 
%  Description: The examples used in reference [1] are:
%
%       Example 7: 3D aircraft stall model for testing stabilization. From
%       Garrard 1977, also studied in Almubarak 2019 ACC paper and others.
%
%       Example 9: Chebychev discretization of Allen-Cahn equation
%
%   Reference: [1] N. A. Corbin and B. Kramer, “Computing solutions to the
%              polynomial-polynomial regulator problem,” in 2024 IEEE 63rd
%              Conference on Decision and Control (CDC), IEEE, Dec. 2024,
%              pp. 2689–2696. doi: 10.1109/cdc56724.2024.10885897.
%
close all; clear; clc;
setKroneckerToolsPath
addpath('examples')

%% Example 7: 3D aircraft stall model
runExample7()

%% Example 9: Allen-Cahn equation
runExample9(33, 4)
% runExample9(129, 4)