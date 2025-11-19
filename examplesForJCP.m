%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Main script to run the examples for JCP %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Author: Nicholas Corbin, UCSD           %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Usage: Run the entire script to generate the examples for the JCP paper 
%         [1]. For the pendulum example (ex. 26), Matlab needs to be able 
%         to call Python to generate the phase portraits. For the 2D 
%         Allen-Cahn example, you will need to install the M-M.E.S.S. 
%         package to access the LR-ADI solver. We also suggest to install
%         the tensor_recursive package for the more efficient tensor
%         Lyapunov solver, or you can use the built-in Bartels-Stewart
%         solver. 
% 
%  Description: The examples used in reference [1] are:
%
%      Example 10: 1D ODE toy example with transcritical bifurcation. The
%      code here prints out the values printed in the paper in Section
%      5.1, namely the open- and closed-loop dynamics for the LQR and PPR
%      controllers. The phase portraits in Fig. 2 in the paper are plotted 
%      using tikz and pgfplots.
%
%      Example 26: Nondimensionalized 2D inverted pendulum. The code here
%      plots the phase portraits, value function contour plots, and HJB
%      residuals for Figs. 3 & 4. A separate function computes the 
%      integrated costs for Table 1.
%
%      Example 9: Chebychev discretization of Allen-Cahn equation with
%      Neumann boundary conditions. The code here compares LQR, SDRE, PPR,
%      and TT-HJB, generating a cost convergence plot Fig. 6 to show if 
%      the controllers are successfully stabilizing the origin, table of 
%      controller costs Table 2, and closed-loop simulation plots Fig. 7. 
%
%      Example 29: Finite element discretization of Allen-Cahn equation 
%      with Neumann boundary conditions. The code computes a reduced-PPR
%      controller and simulates the closed-loop system to demonstrate the
%      scalability of the approach, generating the plots for Fig. 9.
%
%   Reference: [1] N. A. Corbin and B. Kramer, "Nonlinear Feedback Control 
%              in High Dimensions using the Polynomial-Polynomial 
%              Regulator,” submitted.
%
close all; clear; clc;
setKroneckerToolsPath
addpath('examples')
addpath('utils')

%% Example 10: 1D Toy Example with Transcritical Bifurcation
% Prints out the closed-loop dynamics of LQR and PPR controllers
runExample10();

%% Example 26: 2D Inverted Pendulum
% Produces the value function plots, residual plots, and closed-loop phase
% portraits (by calling external python script) for Figs. 3 & 4.
fprintf('\n\nGenerating the plots for Fig. 3 & 4 from the paper... \n\n\n')
for nFterms = 1:2:9
    runExample26(nFterms+1);
end

% Prints out the table of control costs Table 1
runExample26_controlCosts();

%% Example 9: 1D Allen-Cahn equation, Neumann Boundary Conditions
% Produces the cost convergence plot Fig. 6, table of controller costs 
% Table 2, and closed-loop simulation surface plots Fig. 7.
runExample9_Neumann_costConvergencePlot();
runExample9_Neumann();
 
%% Example 29: 2D Allen-Cahn equation, Neumann Boundary Conditions
% Produces the closed-loop simulation plots in Fig. 9.
runExample29(64, 10) % Quick version that runs in about a minute (n=4225)
% runExample29(320, 10) % Paper version (n=103041)
