% function runExample19()
%runExample18 Runs 2D linear example to compare Rekasius result.
%
%   Usage:  runExample19()
%
%   runExample19() runs the default case of a quadratic model from [1]. 
%
%
%   Reference: [1] 
%
%   Part of the PPR repository.
%%

clear; clc; close all; 
set(groot,'defaultLineLineWidth',1.5)

fprintf('Running Example 18\n')

A = [0 1; 0 -1]; B = [0;1];
n=2;
FofXU = @(x,u) (A*x + B*u); 

q4 = zeros(n^4,1); q4(end) = 2; % x2^4 entry; 2 because of 1/2 factor
Q1 = {0,[2 0; 0 0],zeros(n^3,1),q4}; R = 2; 


[ValueFun1, GainsPPR1] = ppr(A, B, Q1, R, 4);
uLin = @(x) (kronPolyEval(GainsPPR1, x, 1));
uCub1 = @(x) (kronPolyEval(GainsPPR1, x));

y = sym('x',[2,1]);
vpa(uCub1(y),3)

%% Plot results to match Bass & Webber Fig. 1
figure('Position', [149 70.3333 786 490.0000]); hold on;

% Linear feedback (which is no feedback)
[t1, X1] = ode45(@(t, x) FofXU(x, uLin(x)), [0, 20], [6.5;0]);
plot(X1(:,1),X1(:,2),'--'); 

% Nonlinear feedback with Φ1, penalizing x2
[t2, X2] = ode45(@(t, x) FofXU(x, uCub1(x)), [0, 20], [6.5;0]);
plot(X2(:,1),X2(:,2)); 

[t2, X2] = ode45(@(t, x) FofXU(x, uCub1(x)), [0, 20], [5.5;0]);
plot(X2(:,1),X2(:,2)); 

[t2, X2] = ode45(@(t, x) FofXU(x, uCub1(x)), [0, 20], [4.5;0]);
plot(X2(:,1),X2(:,2)); 

[t2, X2] = ode45(@(t, x) FofXU(x, uCub1(x)), [0, 20], [3.5;0]);
plot(X2(:,1),X2(:,2)); 

[t2, X2] = ode45(@(t, x) FofXU(x, uCub1(x)), [0, 20], [2.5;0]);
plot(X2(:,1),X2(:,2)); 

[t2, X2] = ode45(@(t, x) FofXU(x, uCub1(x)), [0, 20], [1.5;0]);
plot(X2(:,1),X2(:,2)); 

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
xlim([-1 7]); ylim([-2.8 2]);
% end
