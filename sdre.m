function u = sdre(Ax,Bx,Q,R,x)
%sdre Compute suboptimal control using State-Dependent Riccati Equation.
%
%   Usage: u = sdre(@(x)(A(x)),@(x)(B(x)),Q,R,x)
%
%   Inputs:
%       Ax     - the drift function A(x); can be an anonymous function
%       Bx     - the input map B(x); can be an anonymous function
%       Q      - the state penalty Q
%       R      - the control penalty R
%       x      - the current state value
% 
%   Output:
%       u      - the control u=-Kx, where K is the computed gain matrix
% 
%   Background: For control-affine dynamics
%                   xdot = f(x) + g(x) u 
%               we wish to compute an approximation to the optimal control.
%               Using the State-Dependent Riccati Equation (SDRE) approach,
%               the dynamics are factorized in "semi-linear form" 
%                   xdot = A(x)x + B(x)u 
%               and at each state x, the Riccati equation is solved. 
%
%   Author: Nick Corbin, UCSD
%
%   License: MIT
%
%   Reference: [1] 
%
%  Part of the PPR repository.
%%

% State-dependent system matrices A(x) and B(x)
A = Ax(x); B = Bx(x); 

% Solve the state-dependent Riccati equation: A(x)'P + P*A(x) - P*B(x)*inv(R)*B(x)'P + Q = 0
[P, K, ~, INFO] = icare(A, B, Q, R);

% Check if Riccati equation was solved successfully
if INFO.Report ~= 0
    error('SDRE Riccati equation did not converge.');
end

u = -K*x;
end