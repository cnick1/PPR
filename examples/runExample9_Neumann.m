function runExample9_Neumann()
%runExample9_Neumann Runs the Allen-Cahn example with Nuemman BCs for [1].
%
%   Usage:  runExample9_Neumann()
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
%       In this script, we compute controllers and simulate the closed-loop
%       systems to compare LQR, SDRE, PPR, and TT-HJB, both plotting the
%       solutions and printing the control costs.%
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
uTTHJB = @(z)controlfun_leg(z,-av,av,core2cell(V),gfun,R,umax);
fprintf("completed in %2.2f seconds. \n", toc)
close % close convergence plot
% uTTHJB =  @(z) uLQR(z.'); % For debugging my stuff without running tthjb for so long

%% Simulate closed-loop systems
% Construct original system dynamics
FofXU = @(w,u) (f{1}*w - w.^3 + B*u); % -v.^3 is equivalent to f{3}*x^(3), just faster for simulation purposes

t = [0,1e7];

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-20);

wOffsets = [0.5, 1, 1.5, 2];
for j=1:length(wOffsets)
    clear UxPPR UxPPR_tuned UxPPR_tuned_reduced UxPPR_reduced UxTTHJB UxSDRE
    w0 = wOffsets(j) + cos(2*pi*y).*cos(pi*y); % Modified initial condition

    % Simulate using ode15s and compute performance indexes
    fprintf(" Simulating uncontrolled solution ... "); tic
    [tUnc, Xunc] = ode15s(@(t, v) FofXU(v,uUnc(v)),t, w0, opts);
    fprintf("completed in %2.2f seconds. \n", toc)
    costUnc(j) = trapz(tUnc, sum((Xunc.^2).*diag(Q).', 2));

    fprintf(" Simulating LQR solution ... "); tic
    [tLQR, XLQR] = ode15s(@(t, v) FofXU(v,uLQR(v)),t, w0, opts);
    fprintf("completed in %2.2f seconds. \n", toc)
    UxLQR = uLQR(XLQR.').';
    costLQR(j) = trapz(tLQR, sum((XLQR.^2).*diag(Q).', 2) + R*UxLQR.^2);

    fprintf(" Simulating PPR solution ... "); tic
    [tPPR, XPPR] = ode15s(@(t, v) FofXU(v,uPPR(v)),t, w0, opts);
    fprintf("completed in %2.2f seconds. \n", toc)
    for i=1:length(tPPR); UxPPR(i,1) = uPPR(XPPR(i,:).'); end
    costPPR(j) = trapz(tPPR, sum((XPPR.^2).*diag(Q).', 2) + R*UxPPR.^2);

    fprintf(" Simulating PPR tuned solution ... "); tic
    [tPPR_tuned, XPPR_tuned] = ode15s(@(t, v) FofXU(v,uPPR_tuned(v)),t, w0, opts);
    fprintf("completed in %2.2f seconds. \n", toc)
    for i=1:length(tPPR_tuned); UxPPR_tuned(i,1) = uPPR_tuned(XPPR_tuned(i,:).'); end
    costPPR_tuned(j) = trapz(tPPR_tuned, sum((XPPR_tuned.^2).*diag(Q).', 2) + R*UxPPR_tuned.^2);

    fprintf(" Simulating PPR tuned reduced solution ... "); tic
    [tPPR_tuned_reduced, XPPR_tuned_reduced] = ode15s(@(t, v) FofXU(v,uPPR_tuned_reduced(v)),t, w0, opts);
    fprintf("completed in %2.2f seconds. \n", toc)
    for i=1:length(tPPR_tuned_reduced); UxPPR_tuned_reduced(i,1) = uPPR_tuned_reduced(XPPR_tuned_reduced(i,:).'); end
    costPPR_tuned_reduced(j) = trapz(tPPR_tuned_reduced, sum((XPPR_tuned_reduced.^2).*diag(Q).', 2) + R*UxPPR_tuned_reduced.^2);

    fprintf(" Simulating PPR reduced solution ... "); tic
    [tPPR_reduced, XPPR_reduced] = ode15s(@(t, v) FofXU(v,uPPR_reduced(v)),t, w0, opts);
    fprintf("completed in %2.2f seconds. \n", toc)
    for i=1:length(tPPR_reduced); UxPPR_reduced(i,1) = uPPR_reduced(XPPR_reduced(i,:).'); end
    costPPR_reduced(j) = trapz(tPPR_reduced, sum((XPPR_reduced.^2).*diag(Q).', 2) + R*UxPPR_reduced.^2);

    try
        fprintf(" Simulating SDRE solution ... "); tic
        [tSDRE, XSDRE] = ode15s(@(t, v) FofXU(v,uSDRE(v)),t, w0, opts);
        fprintf("completed in %2.2f seconds. \n", toc)
        for i=1:length(tSDRE); UxSDRE(i,1) = uSDRE(XSDRE(i,:).'); end
        costSDRE(j) = trapz(tSDRE, sum((XSDRE.^2).*diag(Q).', 2) + R*UxSDRE.^2);
    catch
        costSDRE(j) = inf;
    end

    fprintf(" Simulating TT-HJB solution ... "); tic
    [tTTHJB, XTTHJB] = ode15s(@(t, v) FofXU(v,uTTHJB(v.')),t, w0, opts);
    fprintf("completed in %2.2f seconds. \n", toc)
    for i=1:length(tTTHJB); UxTTHJB(i,1) = uTTHJB(XTTHJB(i,:)); end
    costTTHJB(j) = trapz(tTTHJB, sum((XTTHJB.^2).*diag(Q).', 2) + R*UxTTHJB.^2);

    if j == length(wOffsets)
        fprintf('\n# Table 2 Data (Allen-Cahn, Neumann BCs)\n');
        fprintf('# Control costs for different initial condition offsets\n');
        fprintf("      Controller    &     v0=%2.2f      &     v0=%2.2f      &     v0=%2.2f      &     v0=%2.2f     \n",wOffsets)
        fprintf("     %s   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n", 'Uncontrolled', costUnc)
        fprintf("     %s   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n", '     LQR    ', costLQR)
        fprintf("     %s   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n", '     SDRE   ', costSDRE)
        fprintf("     %s   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n", '     PPR    ', costPPR)
        fprintf("     %s   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n", ' PPR reduced', costPPR_reduced)
        fprintf("     %s   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n", ' PPR tuned  ', costPPR_tuned)
        fprintf(" %s  &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n", 'PPR tuned reduced', costPPR_tuned_reduced)
        fprintf("     %s   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n\n", '     TTHJB  ', costTTHJB)
    end
end

% figure
% title('Closed-loop animation')
% for j2 = 1:length(tUnc)
%     tt = tUnc(j2);
%     hold off
%     [~,i] = min(abs(tUnc-tt));
%     plot(Xunc(i,:));
%     hold on
%     [~,i] = min(abs(tLQR-tt));
%     plot(XLQR(i,:));
%     [~,i] = min(abs(tPPR-tt));
%     plot(XPPR(i,:));
%     [~,i] = min(abs(tPPR_tuned-tt));
%     plot(XPPR_tuned(i,:));
%     [~,i] = min(abs(tPPR_reduced-tt));
%     plot(XPPR_reduced(i,:));
%     [~,i] = min(abs(tPPR_tuned_reduced-tt));
%     plot(XPPR_tuned_reduced(i,:));
%     [~,i] = min(abs(tSDRE-tt));
%     plot(XSDRE(i,:)); ylim([-1 1]);
%     [~,i] = min(abs(tTTHJB-tt));
%     plot(XTTHJB(i,:)); ylim([-1 1]);
%     ylim([-3 3]); xlim([0 14])
%     legend('Uncontrolled','LQR','PPR','PPR tuned','PPR reduced','PPR tuned reduced','SDRE','TTHJB')
%     pause(0.005)
% end


set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
tmax = 5; dt = 0.001; t = 0:dt:tmax; % Specify time vector for plotting


% wOffsets = [0.5, 1, 1.5, 2]; % The other cases for the curious
wOffsets = 2;
for j=1:length(wOffsets)
    clear UxPPR UxPPR_tuned UxPPR_tuned_reduced UxPPR_reduced UxTTHJB UxSDRE
    w0 = wOffsets(j) + cos(2*pi*y).*cos(pi*y); % Modified initial condition

    % Simulate using ode15s and compute performance indexes
    fprintf(" Simulating uncontrolled solution ... "); tic
    [tUnc, Xunc] = ode15s(@(t, v) FofXU(v,uUnc(v)),t, w0, opts);
    fprintf("completed in %2.2f seconds. \n", toc)
    costUnc(j) = trapz(tUnc, sum((Xunc.^2).*diag(Q).', 2));

    fprintf(" Simulating LQR solution ... "); tic
    [tLQR, XLQR] = ode15s(@(t, v) FofXU(v,uLQR(v)),t, w0, opts);
    fprintf("completed in %2.2f seconds. \n", toc)
    UxLQR = uLQR(XLQR.').';
    costLQR(j) = trapz(tLQR, sum((XLQR.^2).*diag(Q).', 2) + R*UxLQR.^2);

    fprintf(" Simulating PPR solution ... "); tic
    [tPPR, XPPR] = ode15s(@(t, v) FofXU(v,uPPR(v)),t, w0, opts);
    fprintf("completed in %2.2f seconds. \n", toc)
    for i=1:length(tPPR); UxPPR(i,1) = uPPR(XPPR(i,:).'); end
    costPPR(j) = trapz(tPPR, sum((XPPR.^2).*diag(Q).', 2) + R*UxPPR.^2);

    fprintf(" Simulating PPR tuned solution ... "); tic
    [tPPR_tuned, XPPR_tuned] = ode15s(@(t, v) FofXU(v,uPPR_tuned(v)),t, w0, opts);
    fprintf("completed in %2.2f seconds. \n", toc)
    for i=1:length(tPPR_tuned); UxPPR_tuned(i,1) = uPPR_tuned(XPPR_tuned(i,:).'); end
    costPPR_tuned(j) = trapz(tPPR_tuned, sum((XPPR_tuned.^2).*diag(Q).', 2) + R*UxPPR_tuned.^2);

    fprintf(" Simulating PPR tuned reduced solution ... "); tic
    [tPPR_tuned_reduced, XPPR_tuned_reduced] = ode15s(@(t, v) FofXU(v,uPPR_tuned_reduced(v)),t, w0, opts);
    fprintf("completed in %2.2f seconds. \n", toc)
    for i=1:length(tPPR_tuned_reduced); UxPPR_tuned_reduced(i,1) = uPPR_tuned_reduced(XPPR_tuned_reduced(i,:).'); end
    costPPR_tuned_reduced(j) = trapz(tPPR_tuned_reduced, sum((XPPR_tuned_reduced.^2).*diag(Q).', 2) + R*UxPPR_tuned_reduced.^2);

    fprintf(" Simulating PPR reduced solution ... "); tic
    [tPPR_reduced, XPPR_reduced] = ode15s(@(t, v) FofXU(v,uPPR_reduced(v)),t, w0, opts);
    fprintf("completed in %2.2f seconds. \n", toc)
    for i=1:length(tPPR_reduced); UxPPR_reduced(i,1) = uPPR_reduced(XPPR_reduced(i,:).'); end
    costPPR_reduced(j) = trapz(tPPR_reduced, sum((XPPR_reduced.^2).*diag(Q).', 2) + R*UxPPR_reduced.^2);

    try
        fprintf(" Simulating SDRE solution ... "); tic
        [tSDRE, XSDRE] = ode15s(@(t, v) FofXU(v,uSDRE(v)),t, w0, opts);
        fprintf("completed in %2.2f seconds. \n", toc)
        for i=1:length(tSDRE); UxSDRE(i,1) = uSDRE(XSDRE(i,:).'); end
        costSDRE(j) = trapz(tSDRE, sum((XSDRE.^2).*diag(Q).', 2) + R*UxSDRE.^2);
    catch
        costSDRE(j) = inf;
    end

    fprintf(" Simulating TT-HJB solution ... "); tic
    [tTTHJB, XTTHJB] = ode15s(@(t, v) FofXU(v,uTTHJB(v.')),t, w0, opts);
    fprintf("completed in %2.2f seconds. \n", toc)
    for i=1:length(tTTHJB); UxTTHJB(i,1) = uTTHJB(XTTHJB(i,:)); end
    costTTHJB(j) = trapz(tTTHJB, sum((XTTHJB.^2).*diag(Q).', 2) + R*UxTTHJB.^2);

    
    %%
    fprintf('\n\nPlotting Figure 7 from the paper... \n\n')

    figure('Position',[97 164.3333 1558 577.3333])
    plotIndices = round(linspace(0,1,51).^3*5000+1);
    plotT = t(1:100:end); nplots = length(plotT);
    xx = -1:.025:1; plotdata = [zeros(nplots,length(xx))];
    X = Xunc(plotIndices,:);
    for i=1:length(plotT)
        plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
    end
    subplot(2,4,1)
    mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -0.5 3.05]),
    view(145,35), colormap([0 0 0]); xlabel $z$, ylabel $t$, zlabel $w$
    title("Open-loop"); drawnow

    X = XLQR(plotIndices,:);
    for i=1:length(plotT)
        plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
    end
    subplot(2,4,2)
    mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -0.5 3.05]),
    view(145,35), colormap([0 0 0]); xlabel $z$, ylabel $t$, zlabel $w$
    title("LQR"); drawnow

    X = XSDRE(plotIndices,:);
    for i=1:length(plotT)
        plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
    end
    subplot(2,4,3)
    mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -0.5 3.05]),
    view(145,35), colormap([0 0 0]); xlabel $z$, ylabel $t$, zlabel $w$
    title("SDRE"); drawnow

    X = XTTHJB(plotIndices,:);
    for i=1:length(plotT)
        plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
    end
    subplot(2,4,4)
    mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -0.5 3.05]),
    view(145,35), colormap([0 0 0]); xlabel $z$, ylabel $t$, zlabel $w$
    title("TT-HJB"); drawnow


    X = XPPR(plotIndices,:);
    for i=1:length(plotT)
        plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
    end
    subplot(2,4,5)
    mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -0.5 3.05]),
    view(145,35), colormap([0 0 0]); xlabel $z$, ylabel $t$, zlabel $w$
    title("Degree 5 PPR"); drawnow

    X = XPPR_reduced(plotIndices,:);
    for i=1:length(plotT)
        plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
    end
    subplot(2,4,6)
    mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -0.5 3.05]),
    view(145,35), colormap([0 0 0]); xlabel $z$, ylabel $t$, zlabel $w$
    title("Degree 5 PPR, reduced"); drawnow

    X = XPPR_tuned(plotIndices,:);
    for i=1:length(plotT)
        plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
    end
    subplot(2,4,7)
    mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -0.5 3.05]),
    view(145,35), colormap([0 0 0]); xlabel $z$, ylabel $t$, zlabel $w$
    title("Degree 3 PPR, tuned"); drawnow

    X = XPPR_tuned_reduced(plotIndices,:);
    for i=1:length(plotT)
        plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
    end
    subplot(2,4,8)
    mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -0.5 3.05]),
    view(145,35), colormap([0 0 0]); xlabel $z$, ylabel $t$, zlabel $w$
    title("Degree 3 PPR, tuned \& reduced"); drawnow

    exportgraphics(gcf,sprintf('plots/example9_neumann_v0%2.1f.pdf',wOffsets(j)), 'ContentType', 'vector')

    %%
    fprintf('\n')
end
end

