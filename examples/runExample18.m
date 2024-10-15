% function runExample18()
%runExample18 Runs 3D linear example to check Bass & Webber result.
%
%   Usage:  runExample18()
%
%   runExample18() runs the default case of a quadratic model from [1]. The
%   polynomial version corresponds to the example used in [2].
%
%
%   Reference: [1] 
%
%   Part of the PPR repository.
%%

clear; clc; close all; 
set(groot,'defaultLineLineWidth',1.5)

fprintf('Running Example 18\n')

[A, B, ~] = getSystem18(); n = 3;
FofXU = @(x,u) (A*x + B*u); 

q4 = zeros(n^4,1); q4(41) = 2; % x2^4 entry; 2 because of 1/2 factor
Q1 = {0,zeros(n^2,1),zeros(n^3,1),q4}; R = 1; 

q4 = zeros(n^4,1); q4(81) = 2; % x3^4 entry; 2 because of 1/2 factor
Q2 = {0,zeros(n^2,1),zeros(n^3,1),q4}; 

[ValueFun1, GainsPPR1] = ppr(A, B, Q1, R, 6);
uLin = @(x) (kronPolyEval(GainsPPR1, x, 1));
uCub1 = @(x) (kronPolyEval(GainsPPR1, x));

[ValueFun2, GainsPPR2] = ppr(A, B, Q2, R, 4);
uCub2 = @(x) (kronPolyEval(GainsPPR2, x));

x0 = [5;0;0]; 

[t1, X1] = ode45(@(t, x) FofXU(x, uLin(x)), [0, 15], x0);
[t2, X2] = ode45(@(t, x) FofXU(x, uLin(x)), [0, 15], 4*x0);

[t3, X3] = ode45(@(t, x) FofXU(x, uCub1(x)), [0, 15], x0);
[t4, X4] = ode45(@(t, x) FofXU(x, uCub1(x)), [0, 15], 4*x0);

[t5, X5] = ode45(@(t, x) FofXU(x, uCub2(x)), [0, 15], x0);
[t6, X6] = ode45(@(t, x) FofXU(x, uCub2(x)), [0, 15], 4*x0);


%% Plot results to match Bass & Webber Fig. 1
figure('Position', [1.3143e+03 -44.3333 406.6667 860]); hold on;

% Linear feedback (which is no feedback)
plot(X1(:,2),X1(:,3),':'); plot(X2(:,2),X2(:,3),':')

% Nonlinear feedback with Φ1, penalizing x2
plot(X3(:,2),X3(:,3)); plot(X4(:,2),X4(:,3))

% Nonlinear feedback with Φ2, penalizing x3
plot(X5(:,2),X5(:,3),'--'); plot(X6(:,2),X6(:,3),'--')

xlim([-10 0]); ylim([-16 6]);
% end
