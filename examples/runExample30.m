function runExample30()
%runExample30 Runs the 1D unstable heat equation FEM control example.
%
%   Usage:  runExample30()
%
%   Description: The model, inspired by [1], describes an unstable heat
%   equation problem. The model takes the form
%
%     uₜ(x,t)  = ε uₓₓ(x,t) + λ u(x,t) + μ u(x,t)³
%     uₓ(0,t) = 0                     (Neumann BC)
%     uₓ(1,t) = u(t)                  (Neumann boundary control input)
%
%   where one end is insulated and one end is subject to Neumann boundary
%   control. The FEM model can be written (after multiplying by M⁻¹) as
%
%       ẋ = A x + F₃(x⊗x⊗x) + B u
%       y = C x
%
%   for which we can compute a controller u(x) = K(x) using PPR.
%
%   Reference: [1] D. M. Boskovic, M. Krstic, and W. Liu, "Boundary control
%              of an unstable heat equation via measurement of
%              domain-averaged temperature,” IEEE Transactions on Automatic
%              Control, vol. 46, no. 12, pp. 2022–2028, 2001
%              [2] S. A. Ragab and H. E. Fayed, Introduction to finite
%              element analysis for engineers. Taylor & Francis Group, 2017
%
%   Part of the PPR repository.
%%
fprintf('Running Example 30\n')

% Get dynamics
n = 40; r = n;
fprintf(" Forming FEM model, n=%i ... ",n); tic
[E, f, g, h, xg] = getSystem30(n-1, {0.75,0,0}, {-1,0,1});
fprintf("completed in %2.2f seconds. \n", toc)
g{1} = g{1}(:,2);

% Get value function/controller
q = 0.5; R = 1; degree = 4; options.E = E;
fprintf(" Computing ppr() solution, n=%i, r=%i, d=%i ... ",n,r,4); tic
[~, K] = ppr(f, g, q, R, degree, options);
fprintf(" completed in %2.2f seconds. \n", toc)

uLQR = @(x) kronPolyEval(K, x, 1);
uPPR = @(x) kronPolyEval(K, x, degree-1);

%% Compute controllers
% Setting the cost Q=C.'*C for LR-ADI
% options.lrradi = true;
% nc = 10; nds = round(linspace(1,n,nc));
% C = sparse(1:nc,nds,sqrt(0.5),nc,n); q = C.'*C;
% 
% % Get value function/controller
% R = 1; degree = 4;
% options.C = C; options.E = E;
% options.verbose = false; options.reducedDimension = r;
% fprintf(" Computing ppr() solution w/ lrradi, n=%i, r=%i, d=%i ... ",n,r,4); tic
% [~, K] = ppr(f, g, q, R, degree, options);
% fprintf(" completed in %2.2f seconds. \n", toc)
% 
% uLQR = @(x) kronPolyEval(K, x, 1);
% uPPR = @(x) kronPolyEval(K, x, degree-1);

% To compare with back-stepping controller from [1]
% a = 0.95*pi/2; alpha = 3;
% uBS = @(x) -(alpha + a*tan(a))*x(end) - (alpha*a*tan(a) + a^2/cos(a)^2)*trapz(xg,x);

%% Simulate closed-loop system
% Set up ode function, mass matrix, and Jacobian
% (to maximize efficient sparsity usage)
FofXU = @(x,u) kronPolyEval(f,x) + g{1} * u;
opts_openloop = odeset(Mass=E, Jacobian=f{1});
opts_closloop = odeset(Mass=E, Jacobian=f{1}+g{1}*K{1});

% x0 = 1 - 11*xg.^2 + 18* xg.^3 - 8*xg.^4;
% x0 = 1 - 15*xg.^2 + 26* xg.^3 - 12*xg.^4;
% x0 = 1 - 13.5*xg.^2 + 14* xg.^3 - 1.5*xg.^10;
x0 = 1 - 12*xg.^2 + 12* xg.^3 - 1*xg.^12;
% x0 = 1.1 + 0.*xg;
tmax = 2; t = 0:0.02:tmax; x0 = x0.';% specify for plotting

% Run and time simulations
fprintf(" - Simulating open-loop dynamics ... ");       T0 = tic;
[~, XUNC] = ode15s(@(t, x) FofXU(x,   0    ), t, x0, opts_openloop); fprintf("completed in %2.2f seconds. \n", toc(T0))
fprintf(" - Simulating LQR closed-loop dynamics ... "); T0 = tic;
[~, XLQR] = ode15s(@(t, x) FofXU(x, uLQR(x)), t, x0, opts_closloop); fprintf("completed in %2.2f seconds. \n", toc(T0))
fprintf(" - Simulating PPR closed-loop dynamics ... "); T0 = tic;
[~, XPPR] = ode15s(@(t, x) FofXU(x, uPPR(x)), t, x0, opts_closloop); fprintf("completed in %2.2f seconds. \n", toc(T0))



%% Plot solution
figure('Position',[501 19 632 838])
% figure('Position',[474.3333 340.3333 925.3333 300.6667]);
subplot(3,1,1)
mesh(xg,t(1:size(XUNC,1)),XUNC);
grid on, axis([0 1 0 tmax -1.125 1.125]), view(145,35), colormap([0 0 0]);
xlabel x, ylabel t, zlabel u(x,t), title Uncontrolled; drawnow

% figure('Position',[474.3333 340.3333 925.3333 300.6667]);
subplot(3,1,2)
mesh(xg,t(1:size(XLQR,1)),XLQR);
grid on, axis([0 1 0 tmax -1.125 1.125]), view(145,35), colormap([0 0 0]);
xlabel x, ylabel t, zlabel u(x,t), title LQR; drawnow


subplot(3,1,3)
mesh(xg,t,XPPR);
grid on, axis([0 1 0 tmax -1.125 1.125]), view(145,35), colormap([0 0 0]);
xlabel x, ylabel t, zlabel u(x,t), title PPR; drawnow


% figure
% % Animate
% for i=1:10:length(t)
%     plot(xg,X(i,:));
%     xlim([0 1]); ylim([-1 1]); drawnow
% end



end



