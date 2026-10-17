% function runExample31()
%runExample31 Runs the 1D Fokker-Planck control example to compare with [1].
%
%   Usage:  runExample31()
%
%   Description: The controlled Fokker-Planck equation, as in [1], is given by
%
%     ∂ρ/∂t = νΔρ + ∇·(ρ∇G) + Σᵐⱼ₌₁ u_j ∇·(ρ∇α_j) in Ω×(0,∞),
%     0 = (ν∇ρ + ρ∇G)·n on Γ×(0,∞),
%     ρ(x,0) = ρ0(x) in Ω.
%
%   For more details on the model, see [1]. The finite difference model
%   approximation used here is based on code provided by the authors of
%   [1] via email correspondence.
%
%   Upon discretization and the coordinate changes necessary to center the
%   model, we obtain a finite-dimensional bilinear model of the form
%
%           ẋ = A x + B u + G₁ (x ⊗ u)
%
%   for which we wish to solve an optimal control problem minimizing the
%   cost function
%
%           J(x,u) = ½∫ xᵀ Q x + uᵀ R u dt
%
%   The solution is given by the value function V(x) which solves the HJB
%   PDE, and the optimal control is then
%
%           u(x) = -R⁻¹ gᵀ(x) 𝜕V(x)/𝜕x
%
%   The value function approximation in [1] is in principle identical to
%   that computed in the PPR package, i.e. a Taylor approximation computed
%   via Al'brekht's method. The details of their implementation are
%   slightly different, most notably their use of a quadrature-based k-way
%   Lyapunov solver as opposed to the Bartels-Stewart type solvers the PPR
%   package uses.
%
%   However, the optimal control approximation proposed in [1] is slightly
%   different than the one proposed in the PPR package. Due to the
%   nonlinear input map g(x), it is possible to compute a polynomial
%   feedback law with *some* higher-order terms than the value function
%   approximation. For example, a quadratic value function approximation
%   with Riccati solution V₂ defines a linear feedback law in a Taylor
%   sense: u₁(x) = -R⁻¹ Bᵀ V₂ x. A quadratic *Taylor* approximation of the
%   feedback law would require a cubic value function approximation, as
%   some---but not all---of the quadratic terms in the feedback law contain
%   the cubic value function coefficient. However, due to the nonlinear
%   input g(x), some of the quadratic feedback law terms rely only on lower
%   value function coefficients. One can, therefore, construct *A*
%   quadratic polynomial feedback law using only V₂, simply by including
%   the entire polynomial input map:
%
%       u₂(x) = -R⁻¹ gᵀ(x) V₂ x
%
%   This is what [1] computes, and it is NOT a Taylor approximation (nor do
%   they explicitly claim it is). The PPR package can compute these
%   approximations as well; for bilinear systems, it is straightforward to
%   adapt the code to compute one additional K coefficient. However, it is
%   more efficient to compute g(x), vᵈᵀ(x ⊗ ... ⊗ x), and then multiply
%   them together, so we do this here instead.
%
%   In this script, we will compare a) the true PPR Taylor approximant, b)
%   the Breiten et al. polynomial approximant computed with their code, and
%   c) the modified PPR approximant made to match the Breiten et al.
%   solution. We will compare speed of computation, speed of evaluation,
%   and optimality convergence using the optimal control approximation
%   computed in open-loop also in [1] for comparison.
%
%   References: [1] T. Breiten, K. Kunisch, and L. Pfeiffer, “Numerical study of
%               polynomial feedback laws for a bilinear control problem,”
%               Mathematical Control & Related Fields, vol. 8, no. 3, pp.
%               557–582, 2018, doi: 10.3934/mcrf.2018023.
%
%   Part of the PPR repository.
%%
fprintf('Running Example 31\n')

% Get dynamics
[f, g, h, x0] = getSystem31(false, '2');
q = h{1}.'*h{1}; R = 1e-4;
n = length(x0);
x0 = 2.*x0; 

% Set up time vector for simulations
tmax = 20; t = linspace(0,tmax,2001);

%% Optimal control computed as open-loop signal using gradient descent via a
% two-point boundary value/Pontryagin-like first-order optimality approach
fprintf(" - Simulating optimal closed-loop dynamics via gradient descent ... "); T0 = tic;
sigma= 0.05; theta= 0.7; tolerance= 1e-9;
[xOL, uOL] = optimize(zeros(size(t)),f{1},g{2},g{1},h{1},R,x0,t, ...
    odeset('RelTol',1e-13,'MaxOrder',2),sigma,theta,tolerance);
fprintf(" completed in %2.2f seconds. \n", toc(T0))


%% Use the code from [1] to compute the Breiten et al. polynomial controllers
% To make things clean, I've cleaned up their control computation code.
quad_no = 50;
fprintf(" Computing Breiten et al. controller... "); tic
VBreiten = breitenValueFunction(f{1},g{1},g{2},q,R,quad_no,4);
fprintf(" completed in %2.2f seconds. \n", toc)
uBreiten2 = @(x) -R\(g{1}+g{2}*x).'*computeGain(VBreiten(1:2),x).';
uBreiten3 = @(x) -R\(g{1}+g{2}*x).'*computeGain(VBreiten(1:3),x).';
uBreiten4 = @(x) -R\(g{1}+g{2}*x).'*computeGain(VBreiten,x).';

%% Compute true PPR Controllers (only up to degree 3)
% Get value function/controller
fprintf(" Computing ppr() solution, n=%i d=%i ... ",n,4); tic
[v, K] = ppr(f, g, q, R, 4);
fprintf(" completed in %2.2f seconds. \n", toc)
uLQR = @(x) kronPolyEval(K, x, degree=1);
uPPR2 = @(x) kronPolyEval(K, x, degree=2);
uPPR3 = @(x) kronPolyEval(K, x, degree=3);

% %% Compute reduced PPR Controllers (only up to degree 3)
% % Get value function/controller
% fprintf(" Computing reduced ppr() solution, n=%i d=%i ... ",n,4); tic
% options.reducedDimension = 6;
% [vr, Kr] = ppr(f, g, q, R, 10,options);
% fprintf(" completed in %2.2f seconds. \n", toc)
% uPPR2r = @(x) kronPolyEval(Kr, x, degree=2);
% uPPR3r = @(x) kronPolyEval(Kr, x, degree=3);
% uPPR4r = @(x) kronPolyEval(Kr, x, degree=4);
% uPPR5r = @(x) kronPolyEval(Kr, x, degree=5);
% uPPR6r = @(x) kronPolyEval(Kr, x, degree=6);
% uPPR7r = @(x) kronPolyEval(Kr, x, degree=7);
% uPPR8r = @(x) kronPolyEval(Kr, x, degree=8);
% uPPR9r = @(x) kronPolyEval(Kr, x, degree=9);

%% Construct Breiten-type PPR Controllers by adding the g(x) terms
uPPRb2 = @(x) kronPolyEval(K, x, degree=1) - R\(g{2}*x).'* gradTermEval(v,x,2)*2/2;
uPPRb3 = @(x) kronPolyEval(K, x, degree=2) - R\(g{2}*x).'* gradTermEval(v,x,3)*3/2;
uPPRb4 = @(x) kronPolyEval(K, x, degree=3) - R\(g{2}*x).'* gradTermEval(v,x,4)*4/2;

%% Simulate closed-loop system
FofXU = @(x,u) kronPolyEval(f,x) + kronPolyEval(g,x,scenario='G(x)') * u;
options_ode = odeset('RelTol',1e-8);

% Run and time simulations
% Uncontrolled, LQR, and PPR sims
fprintf(" - Simulating open-loop dynamics ... ");       T0 = tic;
[~, XUNC] = ode15s(@(t, x) FofXU(x,   0    ), t, x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))
fprintf(" - Simulating LQR closed-loop dynamics ... "); T0 = tic;
[~, XLQR] = ode15s(@(t, x) FofXU(x, uLQR(x)), t, x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))
fprintf(" - Simulating PPR2 closed-loop dynamics ... "); T0 = tic;
[~, XPPR2] = ode15s(@(t, x) FofXU(x, uPPR2(x)), t, x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))
fprintf(" - Simulating PPR3 closed-loop dynamics ... "); T0 = tic;
[~, XPPR3] = ode15s(@(t, x) FofXU(x, uPPR3(x)), t, x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))

% % Reduced PPR sims
% fprintf(" - Simulating PPR2r closed-loop dynamics ... "); T0 = tic;
% [~, XPPR2r] = ode15s(@(t, x) FofXU(x, uPPR2r(x)), t, x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))
% fprintf(" - Simulating PPR3r closed-loop dynamics ... "); T0 = tic;
% [~, XPPR3r] = ode15s(@(t, x) FofXU(x, uPPR3r(x)), t, x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))
% fprintf(" - Simulating PPR4r closed-loop dynamics ... "); T0 = tic;
% [~, XPPR4r] = ode15s(@(t, x) FofXU(x, uPPR4r(x)), t, x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))
% fprintf(" - Simulating PPR5r closed-loop dynamics ... "); T0 = tic;
% [~, XPPR5r] = ode15s(@(t, x) FofXU(x, uPPR5r(x)), t, x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))
% fprintf(" - Simulating PPR6r closed-loop dynamics ... "); T0 = tic;
% [~, XPPR6r] = ode15s(@(t, x) FofXU(x, uPPR6r(x)), t, x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))
% fprintf(" - Simulating PPR7r closed-loop dynamics ... "); T0 = tic;
% [~, XPPR7r] = ode15s(@(t, x) FofXU(x, uPPR7r(x)), t, x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))
% fprintf(" - Simulating PPR8r closed-loop dynamics ... "); T0 = tic;
% [~, XPPR8r] = ode15s(@(t, x) FofXU(x, uPPR8r(x)), t, x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))
% fprintf(" - Simulating PPR9r closed-loop dynamics ... "); T0 = tic;
% [~, XPPR9r] = ode15s(@(t, x) FofXU(x, uPPR9r(x)), t, x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))

% Breiten-like PPR sims
fprintf(" - Simulating PPRb2 closed-loop dynamics ... "); T0 = tic;
[~, XPPRb2] = ode15s(@(t, x) FofXU(x, uPPRb2(x)), t, x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))
fprintf(" - Simulating PPRb3 closed-loop dynamics ... "); T0 = tic;
[~, XPPRb3] = ode15s(@(t, x) FofXU(x, uPPRb3(x)), t, x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))
fprintf(" - Simulating PPRb4 closed-loop dynamics ... "); T0 = tic;
[~, XPPRb4] = ode15s(@(t, x) FofXU(x, uPPRb4(x)), t, x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))

% Breiten sims
fprintf(" - Simulating Breiten 2 closed-loop dynamics ... "); T0 = tic;
[~, Xb2] = ode15s(@(t, x) FofXU(x, uBreiten2(x)), t, x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))
fprintf(" - Simulating Breiten 3 closed-loop dynamics ... "); T0 = tic;
[~, Xb3] = ode15s(@(t, x) FofXU(x, uBreiten3(x)), t, x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))
fprintf(" - Simulating Breiten 4 closed-loop dynamics ... "); T0 = tic;
[~, Xb4] = ode15s(@(t, x) FofXU(x, uBreiten4(x)), t, x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))


%% Compute control signals and integrated costs
UxLQR = zeros(length(t),1);
UxPPR2 = zeros(length(t),1);
UxPPR3 = zeros(length(t),1);
UxPPRb2 = zeros(length(t),1);
UxPPRb3 = zeros(length(t),1);
UxPPRb4 = zeros(length(t),1);
Uxb2 = zeros(length(t),1);
Uxb3 = zeros(length(t),1);
Uxb4 = zeros(length(t),1);
for i=1:length(t)
    UxLQR(i,1) = uLQR(XLQR(i,:).');
    UxPPR2(i,1) = uPPR2(XPPR2(i,:).');
    UxPPR3(i,1) = uPPR3(XPPR3(i,:).');
    % UxPPR2r(i,1) = uPPR2r(XPPR2r(i,:).');
    % UxPPR3r(i,1) = uPPR3r(XPPR3r(i,:).');
    % UxPPR4r(i,1) = uPPR4r(XPPR4r(i,:).');
    % UxPPR5r(i,1) = uPPR5r(XPPR5r(i,:).');
    % UxPPR6r(i,1) = uPPR6r(XPPR6r(i,:).');
    % UxPPR7r(i,1) = uPPR7r(XPPR7r(i,:).');
    % UxPPR8r(i,1) = uPPR8r(XPPR8r(i,:).');
    % UxPPR9r(i,1) = uPPR9r(XPPR9r(i,:).');
    UxPPRb2(i,1) = uPPRb2(XPPRb2(i,:).');
    UxPPRb3(i,1) = uPPRb3(XPPRb3(i,:).');
    UxPPRb4(i,1) = uPPRb4(XPPRb4(i,:).');
    Uxb2(i,1) = uBreiten2(Xb2(i,:).');
    Uxb3(i,1) = uBreiten3(Xb3(i,:).');
    Uxb4(i,1) = uBreiten4(Xb4(i,:).');
end

% Compute performance index (cost)
runningCost = @(x,u) 1/2 * sum((x*q) .* x, 2) + 1/2 * R * u.^2;
integratedCostUNC = trapz(t, runningCost(XUNC,0));
integratedCostOL = trapz(t, runningCost(xOL,uOL));
integratedCostLQR = trapz(t, runningCost(XLQR,UxLQR));
integratedCostPPR2 = trapz(t, runningCost(XPPR2,UxPPR2));
integratedCostPPR3 = trapz(t, runningCost(XPPR3,UxPPR3));
% integratedCostPPR2r = trapz(t, runningCost(XPPR2r,UxPPR2r));
% integratedCostPPR3r = trapz(t, runningCost(XPPR3r,UxPPR3r));
% integratedCostPPR4r = trapz(t, runningCost(XPPR4r,UxPPR4r));
% integratedCostPPR5r = trapz(t, runningCost(XPPR5r,UxPPR5r));
% integratedCostPPR6r = trapz(t, runningCost(XPPR6r,UxPPR6r));
% integratedCostPPR7r = trapz(t, runningCost(XPPR7r,UxPPR7r));
% integratedCostPPR8r = trapz(t, runningCost(XPPR8r,UxPPR8r));
% integratedCostPPR9r = trapz(t, runningCost(XPPR9r,UxPPR9r));
integratedCostPPRb2 = trapz(t, runningCost(XPPRb2,UxPPRb2));
integratedCostPPRb3 = trapz(t, runningCost(XPPRb3,UxPPRb3));
integratedCostPPRb4 = trapz(t, runningCost(XPPRb4,UxPPRb4));
integratedCostb2 = trapz(t, runningCost(Xb2,Uxb2));
integratedCostb3 = trapz(t, runningCost(Xb3,Uxb3));
integratedCostb4 = trapz(t, runningCost(Xb4,Uxb4));

%% Print control costs table
% The controllers are grouped in sets of three, corresponding roughly to
% the equivalent computational complexity of each (in the sense, the first
% three compute just a Riccati solution, the second three compute also a
% cubic value function coefficient, the last three use also a quartic value
% function coefficient). The comparison here is therefore between a) a true
% Taylor approximant of degree k to the optimal control, and b) a hybrid polynomial
% controller of degree k+1, using the same information as the degree k
% controller, but hypothesizing that including the full g(x) nonlinearity
% will improve things. In other words, instead of having to add all of the
% degree k+1 terms, is adding just some of them that are cheaper to compute
% a happy medium? The answer is: NO! Adding these terms, in each case,
% HURTS performance. The Breiten controllers are worse than the true Taylor
% approximants of the same order, which maybe could be expected, but they
% are also worse than the Taylor approximants of order 1 less, which simply
% omit the extra term the Breiten controller adds!

% Note: Since the PPR "B" controllers are identical to the Breiten controllers,
% we could omit them
fprintf('\n\n# Table Data (Fokker-Planck Control Costs)\n');
fprintf('# Integrated costs for different controllers; the extra terms \n# in the Breiten controllers consistently hurt performance\n');
fprintf('    Controller         &  Integrated Cost\n');
fprintf('    %-18s & %12.6f\n', 'Uncontrolled', integratedCostUNC);
fprintf('    %-18s & %12.6f\n', 'LQR', integratedCostLQR);
% fprintf('    %-18s & %12.6f\n', 'PPRb2', integratedCostPPRb2);
fprintf('    %-18s & %12.6f\n', 'Breiten 2', integratedCostb2);
fprintf('    %-18s & %12.6f\n', 'PPR2', integratedCostPPR2);
% fprintf('    %-18s & %12.6f\n', 'PPR2r', integratedCostPPR2r);
% fprintf('    %-18s & %12.6f\n', 'PPRb3', integratedCostPPRb3);
fprintf('    %-18s & %12.6f\n', 'Breiten 3', integratedCostb3);
fprintf('    %-18s & %12.6f\n', 'PPR3', integratedCostPPR3);
% fprintf('    %-18s & %12.6f\n', 'PPR3r', integratedCostPPR3r);
% fprintf('    %-18s & %12.6f\n', 'PPRb4', integratedCostPPRb4);
fprintf('    %-18s & %12.6f\n', 'Breiten 4', integratedCostb4);
% fprintf('    %-18s & %12.6f\n', 'PPR4r', integratedCostPPR4r);
% fprintf('    %-18s & %12.6f\n', 'PPR5r', integratedCostPPR5r);
% fprintf('    %-18s & %12.6f\n', 'PPR6r', integratedCostPPR6r);
% fprintf('    %-18s & %12.6f\n', 'PPR7r', integratedCostPPR7r);
% fprintf('    %-18s & %12.6f\n', 'PPR8r', integratedCostPPR8r);
% fprintf('    %-18s & %12.6f\n', 'PPR9r', integratedCostPPR9r);
fprintf('    %-18s & %12.6f\n', 'Optimal (open-loop)', integratedCostOL);
fprintf('\n');

% for i=1:length(t)
%     hold off;
%     plot(XUNC(i,:),'DisplayName','Uncontrolled'); hold on;
%     plot(XLQR(i,:),'DisplayName','LQR')
%     plot(xOL(i,:),'DisplayName','Optimal')
%     plot(XPPR3(i,:),':','DisplayName','PPR')
%     plot(XPPRb4(i,:),'r:','DisplayName','Breiten 4')
%     legend
%     pause(0.01);
% end

%%
fprintf('\n\nPlotting Figure from the paper... \n\n')
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
figure('Position',[293.8000 462.6000 754.2000 295.4000]); tiledlayout(1,2);
dt = 25; dx = 2; xx = 1:n;%linspace(-6,6,n);
% Open-loop uncontrolled plot
nexttile;
mesh(xx(1:dx:end),t(1:dt:end),XUNC(1:dt:end,1:dx:end))
grid on, axis([1 n 0 5 -0.5 1.05]),
view(35,35), colormap([0 0 0]); xlabel $i$, ylabel $t$, zlabel $x_i$
title("Uncontrolled"); drawnow
% Open-loop optimal control plot
nexttile;
mesh(xx(1:dx:end),t(1:dt:end),xOL(1:dt:end,1:dx:end))
grid on, axis([1 n 0 5 -0.5 1.05]),
view(35,35), colormap([0 0 0]); xlabel $i$, ylabel $t$, zlabel $x_i$
title("Open-loop Optimal Control"); drawnow
exportgraphics(gcf,'plots/example31_evolution.pdf', 'ContentType', 'vector')

% fprintf('\n\nPlotting Figure from the paper... \n\n')
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% figure('Position',[293.8000 462.6000 754.2000 295.4000]); tiledlayout(1,4);
% dt = 25; dx = 2; xx = linspace(-6,6,n);
% % Open-loop uncontrolled plot
% nexttile;
% mesh(xx(1:dx:end),t(1:dt:end),XUNC(1:dt:end,1:dx:end))
% grid on, axis([-6 6 0 5 -0.5 1.05]),
% view(145,35), colormap([0 0 0]); xlabel $z$, ylabel $t$, zlabel $w$
% title("Uncontrolled"); drawnow
% % Open-loop optimal control plot
% nexttile;
% mesh(xx(1:dx:end),t(1:dt:end),xOL(1:dt:end,1:dx:end))
% grid on, axis([-6 6 0 5 -0.5 1.05]),
% view(145,35), colormap([0 0 0]); xlabel $z$, ylabel $t$, zlabel $w$
% title("Open-loop Optimal Control"); drawnow
% % LQR control plot
% nexttile;
% mesh(xx(1:dx:end),t(1:dt:end),XLQR(1:dt:end,1:dx:end))
% grid on, axis([-6 6 0 5 -0.5 1.05]),
% view(145,35), colormap([0 0 0]); xlabel $z$, ylabel $t$, zlabel $w$
% title("LQR Control"); drawnow
% % Breiten2 control plot
% nexttile;
% mesh(xx(1:dx:end),t(1:dt:end),Xb2(1:dt:end,1:dx:end))
% grid on, axis([-6 6 0 5 -0.5 1.05]),
% view(145,35), colormap([0 0 0]); xlabel $z$, ylabel $t$, zlabel $w$
% title("Breiten et al. Quadratic Control"); drawnow

%% Plot control signals
figure; hold on;
set(0, 'DefaultLineLineWidth', 1.5);
plot(t,Uxb2,'r-','DisplayName','$u_{b2}(t)$') % Breiten quadratic
plot(t,Uxb3,'g','DisplayName','$u_{b3}(t)$') % Breiten cubic
plot(t,Uxb4,'b','DisplayName','$u_{b4}(t)$') % Breiten quartic
plot(t,UxPPRb2,'y--','DisplayName','$u_{bp2}(t)$') % Breiten quadratic computed with PPR package
plot(t,UxPPRb3,'m--','DisplayName','$u_{bp3}(t)$') % Breiten cubic computed with PPR package
plot(t,UxPPRb4,'c--','DisplayName','$u_{bp4}(t)$') % Breiten quartic computed with PPR package
plot(t,UxLQR,'k:','DisplayName','$u_1(t)$') % LQR (linear PPR)
plot(t,UxPPR2,'r:','DisplayName','$u_2(t)$') % Quadratic PPR
plot(t,uOL,'k','DisplayName','$u_{opt}(t)$') % True optimal control
plot(t,UxPPR3,'g:','DisplayName','$u_3(t)$') % Cubic PPR
% plot(t,UxPPR2r,'r:','DisplayName','u_{2r}(t)')
% plot(t,UxPPR3r,'g:','DisplayName','u_{3r}(t)')
% plot(t,UxPPR4r,'b:','DisplayName','u_{4r}(t)')
% plot(t,UxPPR5r,'m:','DisplayName','u_{5r}(t)')
% plot(t,UxPPR6r,'c:','DisplayName','u_{6r}(t)')
% plot(t,UxPPR7r,'y:','DisplayName','u_{7r}(t)')
% plot(t,UxPPR8r,'k:','DisplayName','u_{8r}(t)')
% plot(t,UxPPR9r,'r:','DisplayName','u_{9r}(t)')
legend
xlim([0 5])


figure('Position',[293.8000 462.6000 754.2000 295.4000]); tiledlayout(1,2);
nexttile;
hold on;

plot(t,uOL,'Color','k','LineStyle','-','LineWidth',1.8,...
    'DisplayName','$u_{opt}(t)$');
plot(t,Uxb2,'b','LineWidth',1.2,'DisplayName','$u_{b2}(t)$');
plot(t,UxPPRb2,'--', 'Color', [0.75, 0.75, 0.75],'LineWidth',1,...
    'Marker','o','MarkerIndices',1:40:length(t),'MarkerSize',5,...
    'DisplayName','$u_{b2}^{\mathrm{PPR}}(t)$');
plot(t,Uxb3,'r','LineWidth',1.2,'DisplayName','$u_{b3}(t)$');
plot(t,UxPPRb3,'--', 'Color', [0.75, 0.75, 0.75],'LineWidth',1,...
    'Marker','s','MarkerIndices',1:40:length(t),'MarkerSize',6,...
    'DisplayName','$u_{b3}^{\mathrm{PPR}}(t)$');
plot(t,Uxb4,'g','LineWidth',1.2,'DisplayName','$u_{b4}(t)$');
plot(t,UxPPRb4,'--', 'Color', [0.75, 0.75, 0.75],'LineWidth',1,...
    'Marker','^','MarkerIndices',1:40:length(t),'MarkerSize',5,...
    'DisplayName','$u_{b4}^{\mathrm{PPR}}(t)$');
% plot(t,UxLQR,'b:','LineWidth',1.4,'DisplayName','$u_1(t)$');
% plot(t,UxPPR2,'r--','LineWidth',1.4,'DisplayName','$u_2(t)$');
% plot(t,UxPPR3,'g-.','LineWidth',1.4,'DisplayName','$u_3(t)$');
legend;
xlim([0 5]); ylim([-5 35])
xlabel('$t$'); ylabel('$u(t)$');
title('Breiten et al. Controllers')

nexttile;
hold on;

plot(t,uOL,'Color','k','LineStyle','-','LineWidth',1.8,...
    'DisplayName','$u_{opt}(t)$');
plot(t,UxLQR,'b:','LineWidth',1.4,'DisplayName','$u_1(t)$');
plot(t,UxPPR2,'r--','LineWidth',1.4,'DisplayName','$u_2(t)$');
plot(t,UxPPR3,'g-.','LineWidth',1.4,'DisplayName','$u_3(t)$');

legend;
xlim([0 5]); ylim([-5 35])
xlabel('$t$'); ylabel('$u(t)$');
title('PPR Controllers')

exportgraphics(gcf,'plots/example31_controlSignals.pdf', 'ContentType', 'vector')

% end



function y = gradTermEval(v,x,k)
% Compute the k-1th degree term of the gradient of the value function (the
% derivative term corresponding to the kth degree term in the value
% function). Due to the symmetry of vk, we can do this by reshaping and
% multiplying by x, k-1 times, contracting the tensor Vk down to a vector
% as we go.
n = size(x, 1);
y = v{k};
for i = k:-1:2
    % Using kron-vec identity (level-3 vs level-2 BLAS)
    y = reshape(y,n^(i-1),n) * x;
end
end

%% The code below was provided by Breiten et al.
function gain = computeGain(T,x)
% Note: this is moreso just the gradient of the value function, not really
% the gain matrix. This is their implementation, which is not as efficient
% as what we do (you can recycle xt instead of starting from scratch every
% time, and instead of forming x^k you can just multiply with x k-times
% (which is what gradTermEval above does).
% Their T is akin to our v, but our v's are symmetric vectors (vectorized
% symmetric tensors) whereas it is not clear without digging further what
% properties their T's satisfy, so I will leave their implementation intact.
n = size(x, 1);
gain = 0;
for k = 2:size(T,2)
    xt = x';
    for j = 1:k-2
        xt = kron(x',xt);
    end
    gain = gain + 1/factorial(k-1)*xt*reshape(T{k},n^(k-1),n);
end
end

function T = breitenValueFunction(A,B,G1,q,R,quad_no,degree)
% This is a cleaned up version of the code from [1], packaged as a function
T{2} = icare(full(A),full(B),full(q),R);
Api = A-R\B*B'*T{2};

for j = 3:degree
    rhs  = set_up_rhs(T,B,G1,j,R);
    
    % Approximation by Quadrature
    T{j} = solv_kron_ls( Api',rhs,j,quad_no );
end

end

function [ rhs ] = set_up_rhs( T,B,N,k,gamma )
%SET_UP_RHS Summary of this function goes here
%   Detailed explanation goes here
n = size(B,1);
u{1} = mat_mult( T{2},2,n,1,B' );
for j = 2:k-2
    u{j} = reshape(mat_mult( T{j+1},j+1,n,1,B' ),n*ones(1,j));
    for l = 1:j
        u{j} = u{j} + mat_mult( T{j},j,n,l,N' );
    end
end

u{k-1} = 0;
for l = 1:k-1
    u{k-1} = u{k-1} + mat_mult( T{k-1},k-1,n,l,N' );
end

rhs = 0;
for i = 1:k-1
    v = kron(u{i}(:),u{k-i}(:));
    v = reshape(v,n*ones(1,k));
    
    p = nchoosek(1:k,i);
    clear q
    for l = 1:size(p,1)
        q(l,:) = setdiff(1:k,p(l,:));
    end
    PQ = [q,p];
    for j = 1:size(p,1)
        perv = ipermute(v,PQ(j,:));
        rhs = rhs + perv;
    end
end
rhs = rhs/(2*gamma);

end

function [ x ] = solv_kron_ls( A,f,d,k )
%COMPUTE_INVERSE Summary of this function goes here
%   Detailed explanation goes here

lam_min_A = -max(real(eig(A)));
lam_min = d*lam_min_A;
x = 0;
hst = pi/sqrt(k);
for j = -k:k
    tj = log(exp(j*hst)+sqrt(1+exp(2*j*hst)));
    wj = hst/sqrt(1+exp(-2*j*hst));
    D = expm(2*tj/lam_min*A);
    x = x-2*wj/lam_min*mat_mult_all( f,D,d );
    % j
end


end

function [x, u] = optimize(uInit,A,N,B,C,gamma,xInit,tsim,options,sigma,beta,tolerance)
figure;
u= uInit;
tsimBack= fliplr(tsim);
[~,x] = ode15s(@state_equation,tsim,xInit,options,A,N,B,u,tsim);
x= x';
[~,p] = ode15s(@costate_equation,tsimBack,zeros(length(xInit),1),options,A,N,C,x,u,tsim);
p= p';
p = fliplr(p);
gradU = gamma*u + B'*p + ones(1,size(x,1))*(p.*(N*x));
error= L2norm_gradient(gradU, tsim);
test= (L2norm_gradient(gradU, tsim)< tolerance);

costSaved= costFunction(u,A,N,B,C,gamma,xInit,tsim,options);
k= 0;

while (~test)&&(k<1000)
    k= k+1;
    step_length = compute_step_length(u,gradU,A,N,B,C,gamma,xInit,tsim,options,sigma,beta);
    u= u - step_length*gradU;
    [~,x] = ode15s(@state_equation,tsim,xInit,options,A,N,B,u,tsim);
    x= x';
    [~,p] = ode15s(@costate_equation,tsimBack,zeros(length(xInit),1),options,A,N,C,x,u,tsim);
    p= p';
    p = fliplr(p);
    gradU = gamma*u + B'*p + ones(1,size(x,1))*(p.*(N*x));
    error= L2norm_gradient(gradU, tsim); semilogy(k,error,'*'); hold on; grid on;
    ylabel('Error'); xlabel('Iteration'); title('Gradient Descent Convergence')
    drawnow
    test= (L2norm_gradient(gradU, tsim)< tolerance);
end
close
costNew= costFunction(u,A,N,B,C,gamma,xInit,tsim,options);
costSaved-costNew;

x=x.';
u=u.';

end

function step_length = compute_step_length(u,gradU,A,N,B,C,gamma,xInit,tsim,options,sigma,beta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

currentCost = costFunction(u,A,N,B,C,gamma,xInit,tsim,options);
slope = L2norm_gradient(gradU, tsim);
step_length= 500;
newCost= costFunction(u-step_length*gradU,A,N,B,C,gamma,xInit,tsim,options);
test= (newCost <= currentCost - sigma*step_length*slope);

while (~test)
    step_length= beta*step_length;
    newCost= costFunction(u-step_length*gradU,A,N,B,C,gamma,xInit,tsim,options);
    test= (newCost <= currentCost - sigma*step_length*slope);
end

variation= newCost-currentCost;


end

function dx = state_equation(t,x,A,N,B,uvalues,tsim)

u = interp1(tsim,uvalues',t);
dx = A*x + N*x*u +B*u;
%t
end

function dp = costate_equation(t,p,A,N,C,xvalues,uvalues,tsim)

x = interp1(tsim,xvalues',t)';
u = interp1(tsim,uvalues',t);
dp = -((C'*C)*x+ (A'+(u*N'))*p);
%t
end

function slope = L2norm_gradient(gradU, tsim)

fun= @(t) (interp1(tsim,gradU',t)).^2;
slope = integral(fun,0,tsim(end));

end

function cost = costFunction(u,A,N,B,C,gamma,xInit,tsim,options)

[~,x] = ode23(@state_equation,tsim,xInit,options,A,N,B,u,tsim);
x= x';
fun= @(t) integral_cost(t,C,gamma,x,u,tsim);
cost= integral(fun,0, tsim(end));

end

function [ TA ] = mat_mult( T,d,n,mu,A )
%MAT_MULT Summary of this function goes here
%   Detailed explanation goes here

m = size(A,1);
TA = reshape(T,n*ones(1,d));
TA = A*reshape(permute(TA,[mu,setdiff(1:d,mu)]),n,n^(d-1));
TA = reshape(TA,[m,n*ones(1,d-1)]);
TA = ipermute(TA,[mu,setdiff(1:d,mu)]);

end


function [ TA ] = mat_mult_all( T,A,d )
%MAT_MULT Summary of this function goes here
%   Detailed explanation goes here

m = size(A,1);
n = size(A,2);

TA = T;
for i = 1:d
    TA = A*reshape(TA,n,numel(TA)/n);
    TA = reshape(TA,[m*ones(1,i),n*ones(1,d-i)]);
    TA = permute(TA,[2:d,1]);
end

end

function z = integral_cost(t,C,gamma,xvalues,uvalues,tsim)

u = interp1(tsim,uvalues',t);
x = interp1(tsim,xvalues',t)';

z= 0.5*(ones(1,size(C,1))*((C*x).^2) + gamma*u.^2);

end

