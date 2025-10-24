function runExample10()
%runExample10 Runs the 1D transcritical bifurcation example for [1].
%
%   Usage: runExample10()
%
%   Background: Consider the following 1D ODE model: 
%
%       ẋ = 1/10 x - 2 x² + x³ + u
%       y = x
%
%   We compute the value function and a linear and cubic approximation of
%   the optimal feedback law. The closed-loop dynamics are printed to the
%   command window. 
%
%   Reference: [1] N. A. Corbin and B. Kramer, "Nonlinear Feedback Control 
%                  in High Dimensions using the Polynomial-Polynomial 
%                  Regulator,” in preparation.
%
%   Part of the PPR repository.
%%
fprintf('\nRunning Example 10, 1D transcritical bifurcation example \n')

degree = 4;
[f, g, ~] = getSystem10(); m = 1; n = 1;
Q2 = 1; q = {[],Q2,2/5,0}; R = 1;

% Full PPR solution (LQR is just the first term)
fprintf("Computing ppr() solution, n=%i, d=%i ... ", n, degree); tic
[~, GainsPPR] = ppr(f, g, q, R, degree);
fprintf("completed in %2.2f seconds. \n", toc)

% Construct control laws
uOpenLoop = @(z) zeros(m,1);
uLQR = @(z) (kronPolyEval(GainsPPR, z, 1));
uPPR = @(z) (kronPolyEval(GainsPPR, z));

%% Simulate closed-loop systems
% Construct original system dynamics
F = @(x) kronPolyEval(f, x);
G = @(x) (g{1} + kronPolyEval(g(2:end), x));
FofXU = @(x,u) (F(x) + G(x)*u);

syms x
fprintf('  Open-loop dynamics:   dx/dt = %s\n', char(vpa(FofXU(x, uOpenLoop(x)),4)))
fprintf('        LQR dynamics:   dx/dt = %s\n', char(vpa(FofXU(x, uLQR(x)),     4)))
fprintf('        PPR dynamics:   dx/dt = %s\n', char(vpa(FofXU(x, uPPR(x)),     4)))

end
