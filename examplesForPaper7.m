%  A script to run the examples in PPR Pade/Chisholm/Canterbury paper
%  The examples are the following:
%
%     Example 2: a 2-DOF quadratic-bilinear model based on
%     Scherpen/Kawano [2]
% 
%     Example 11: 2D inverted pendulum. We can choose different degree
%     approximations to the dynamics, different degree energy functions
%     to compute, etc.
% 
%     Example 12: 2D academic/toy model from Fujimoto/Scherpen
% 
%     Example 14: 2D double pendulum gradient system
%
close all; clear; clc;
setKroneckerToolsPath
addpath('examples')

exportData = false;

%% Example 2: 2D QB model
% Ell = 3
runExample2_chisholmApproximant(3,2,1) 
% runExample2_chisholmApproximant(3,3,2)

% runExample2_chisholmApproximant(5,4,1)
% runExample2_chisholmApproximant(5,5,2)

%% Example 12: 2D toy Scherpen model
% Ell = 3
runExample12_chisholmApproximant(3,2,1) 

% Ell = 7
runExample12_chisholmApproximant(7,6,1)
runExample12_chisholmApproximant(7,7,2)

% These seem to not work
% runExample12_chisholmApproximant(7,4,1)
% runExample12_chisholmApproximant(7,5,2)
% runExample12_chisholmApproximant(7,6,3) 
% runExample12_chisholmApproximant(7,7,4)

% runExample12_chisholmApproximant(7,2,1)
% runExample12_chisholmApproximant(7,3,2)
% runExample12_chisholmApproximant(7,4,3)
% runExample12_chisholmApproximant(7,5,4)
% runExample12_chisholmApproximant(7,6,5)
% runExample12_chisholmApproximant(7,7,6)

runExample12_chisholmApproximant(5,3,4)
runExample12_chisholmApproximant(5,4,5)
runExample12_chisholmApproximant(5,5,6)
runExample12_chisholmApproximant(5,6,5)
runExample12_chisholmApproximant(7,3,8)
runExample12_chisholmApproximant(7,4,5)
runExample12_chisholmApproximant(7,4,9)
runExample12_chisholmApproximant(7,6,4)

for i=7:2:11
    for j=2:i+2 
        for k=1:i+2
            try 
                runExample12_chisholmApproximant(i,j,k)
                fprintf('\n                                                                                 Success with %i, %i, %i\n\n',i,j,k)
                pause(3)
                close all
            catch

            end
        end
    end
end


%% Example 14: 2D double pendulum gradient system
% Ell = 3
runExample14_chisholmApproximant(3,2,1) 
% runExample14_chisholmApproximant(3,3,2)

%% Example 11: 2D inverted pendulum
% Ell = 3
runExample11_chisholmApproximant(3,2,1) 
% runExample11_chisholmApproximant(3,3,2)
% 
% % Ell = 5
% runExample11_chisholmApproximant(5,4,1)
% runExample11_chisholmApproximant(5,5,2)
% 
% runExample11_chisholmApproximant(5,2,1)
% runExample11_chisholmApproximant(5,3,2)
% runExample11_chisholmApproximant(5,4,3)
% runExample11_chisholmApproximant(5,5,4) %
% 
% % Ell = 7
% runExample11_chisholmApproximant(7,6,1)
% runExample11_chisholmApproximant(7,7,2)
% 
% runExample11_chisholmApproximant(7,4,1)
% runExample11_chisholmApproximant(7,5,2)
% runExample11_chisholmApproximant(7,6,3) %
% runExample11_chisholmApproximant(7,7,4)
% 
% runExample11_chisholmApproximant(7,2,1)
% runExample11_chisholmApproximant(7,3,2)
% runExample11_chisholmApproximant(7,4,3)
% runExample11_chisholmApproximant(7,5,4)
% runExample11_chisholmApproximant(7,6,5)
% runExample11_chisholmApproximant(7,7,6)




