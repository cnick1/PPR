function runExample9_controlEvaluationTiming()
%runExample9_controlEvaluationTiming Runs the Allen-Cahn example with Nuemman BCs for [1].
%
%   Usage:  runExample9_controlEvaluationTiming()
%
%   Background: The Allen-Cahn equation is
%                   uₜ = eps*uₓₓ+u-u³, uₓ(-1)=0, uₓ(1)=0
%       The spatial domain is discretized using Chebychev points, and the
%       spatial derivative becomes uₓ = D*u with differentiation matrix D.
%       The spatial domain x and differentiation matrix D are given by the
%       cheb() function from [3]. The system has equilibrium solutions
%       u(x) = -1, 0, 1, which all satisfy the Neumann boundary conditions.
%       The family of steady-state solutions given by u(x) = tanh((x-x0)/sqrt(2*eps))
%       also form equilibrium solutions.
%
%       Upon discretization, the model is
%                   uₜ = eps*D²*u+u-u³
%
%           -> A = eps*D² + I
%       The equilibrium at the origin is the unstable u(x)=0, so we seek to
%       stabilize it, as was done in [2].
%       
%       In this script, we compute controllers evaluate the online cost of
%       evaluating the controller many times. 
%
%   References: [1] N. A. Corbin and B. Kramer, "Nonlinear Feedback Control
%               in High Dimensions using the Polynomial-Polynomial
%               Regulator,” in preparation.
%               [2] S. Dolgov, D. Kalise, and K. Kunisch, "Tensor
%               decomposition methods for high-dimensional
%               Hamilton-Jacobi-Bellman equations," SIAM Journal on
%               Scientific Computing, vol. 43, no. 3, pp. A1625–A1650, Jan.
%               2021, doi: 10.1137/19m1305136.
%               [3] L. N. Trefethen, Spectral methods in MATLAB. Society
%               for Industrial and Applied Mathematics, 2000.
%               doi: 10.1137/1.9780898719598.
%
%   Part of the PPR repository.
%%
n = 14;
eps = 0.5;
r = 3;

fprintf('\nRunning Example 9, Allen-Cahn example with Nuemman BCs, for ε=%.1f\n',eps)

%% Get dynamics and define control problem
[f, B, Q, ~, y] = getSystem9Neumann(eps, n);
R = 1e-1; m = size(B,2);

%% Construct controllers
% Open-loop (uncontrolled) controller
uUnc = @(z) zeros(m,1);

% LQR Controller
fprintf(" Computing lqr solution, n=%i, d=%i ... ",n,2); tic
options = struct; options.verbose = false;
[~, GainsPPR] = ppr(f, B, Q, R, 2, options);
fprintf("completed in %2.2f seconds. \n", toc)

% PPR Controller
fprintf(" Computing ppr() solution, n=%i, d=%i ... ",n,6); tic
options = struct; options.verbose = false;
[~, GainsPPR] = ppr(f, B, Q, R, 6, options);
fprintf("completed in %2.2f seconds. \n", toc)

uPPR = @(z) (kronPolyEval(GainsPPR, z));

% Tuned PPR Controller
fprintf(" Computing tuned ppr() solution, n=%i, d=%i ... ",n,4); tic
options = struct; options.verbose = false; % It appears that the optimal control strategy for larger offset initial conditions is actually to temper the controller; use a little bit of input but sort of ride it out for a while before kicking in further. Can I use that as intuition to choose a higher order Q4 or R *artificially* to solve the quadratic cost problem with some insight?
[~, GainsPPR_tuned] = ppr(f, B, {0,Q,0,0.25}, {R,0.009,0}, 4, options);
fprintf("completed in %2.2f seconds. \n", toc)

uPPR_tuned = @(z) (kronPolyEval(GainsPPR_tuned, z));

% Reduced PPR Controller
fprintf(" Computing ppr() solution, n=%i, r=%i, d=%i ... ",n,r,6); tic
options = struct; options.verbose = false; options.reducedDimension = r; options.h = B.';
[~, GainsPPR_reduced, options] = ppr(f, B, Q, R, 6, options);
fprintf("completed in %2.2f seconds. \n", toc)

uPPR_reduced = @(z) (kronPolyEval(GainsPPR_reduced, z));

% Tuned Reduced PPR Controller
fprintf(" Computing tuned reduced ppr() solution, n=%i, r=%i, d=%i ... ",n,r,4); tic
options = struct; options.verbose = false; options.reducedDimension = r; options.h = B.';
[~, GainsPPR_tuned_reduced, options] = ppr(f, B, {0,Q,0,0.25}, {R,0.009,0}, 4, options);
fprintf("completed in %2.2f seconds. \n", toc)

uPPR_tuned_reduced = @(z) (kronPolyEval(GainsPPR_tuned_reduced, z));

% LQR Controller (first term in PPR controller)
uLQR = @(z) (kronPolyEval(GainsPPR, z, 1));

% SDRE Controller
uSDRE = @(z) sdre(@(y)(f{1} - diag(y.^2)),@(y)(B),Q,R,z);

% TT-HJB Controller
addpath(genpath('../TT-Toolbox/'),genpath('../tamen/'),genpath('../TT-HJB/'))
rng("default")
umax = inf;nv = 5;av = 3;tol = 1e-3;mu = 100;
% Prepare handles for ODE component functions for HJB
ffun = @(x,i)x*f{1}(i,:).' - x(:,i).^3; % nonlinear Dynamics
gfun = @(x,i)B(i)*ones(size(x,1),1); % Actuator (here just const vector)
lfun = @(x)sum((x.^2).*diag(Q).', 2);
% Call Dolgov code
fprintf(" Computing TT-HJB solution ... "); figure; tic
[printout, V] = evalc('hjb_leg(size(f{1},1), nv, av, ffun, gfun, lfun, R, tol, mu, [], umax)'); % Call with evalc to keep command window clean
uTTHJB = @(z)controlfun_leg(z.',-av,av,core2cell(V),gfun,R,umax);
fprintf("completed in %2.2f seconds. \n", toc)
close % close convergence plot
% uTTHJB =  @(z) uLQR(z); % For debugging my stuff without running tthjb for so long

%% Evaluate many times 
w0 = 1 + cos(2*pi*y).*cos(pi*y); % Modified initial condition
N = 10000;
tic
for i=1:N
 u = uUnc(w0);
end
tUnc = toc;

tic
for i=1:N
 u = uLQR(w0);
end
tLQR = toc;

tic
for i=1:N
 u = uPPR(w0);
end
tPPR = toc;

tic
for i=1:N
 u = uTTHJB(w0);
end
tTTHJB = toc;

tic
for i=1:N
 u = uSDRE(w0);
end
tSDRE = toc;

tic
for i=1:N
 u = uPPR_reduced(w0);
end
tPPR_reduced = toc;

tic
for i=1:N
 u = uPPR_tuned(w0);
end
tPPR_tuned = toc;

tic
for i=1:N
 u = uPPR_tuned_reduced(w0);
end
tPPR_tuned_reduced = toc;




%% 
fprintf('\n# Table 2 Data (Allen-Cahn, Neumann BCs)\n');
fprintf('# Time to evaluate the control law 10,000 times. \n');
fprintf("      Controller        &         CPU-Time     \n")
fprintf("     %s       &  %13.3f   \n", 'Uncontrolled', tUnc)
fprintf("     %s       &  %13.3f   \n", '    LQR     ', tLQR)
fprintf("     %s       &  %13.3f   \n", '   SDRE     ', tSDRE)
fprintf("    %s   &  %13.3f   \n", 'PPR (degree 6)   ', tPPR)
fprintf("     %s       &  %13.3f   \n", 'PPR reduced ', tPPR_reduced)
fprintf(" %s  &  %13.3f   \n", ' PPR tuned (degree 4)', tPPR_tuned)
fprintf(" %s      &  %13.3f   \n", 'PPR tuned reduced', tPPR_tuned_reduced)
fprintf("     %s       &  %13.3f   \n\n", '     TTHJB  ', tTTHJB)
end 