function runExample9_Neumann(n, degree, eps, r)
%runExample9_Neumann Runs the Allen-Cahn example with Nuemman BCs.
%
%   Usage:  runExample9_Neumann(n,degree,eps,r)
%
%   Inputs: n       - desired state dimension
%           degree  - desired polynomial degree of value function to compute
%           eps     - diffusion coefficient
%           r       - ROM dimension; if r=n, no MOR is performed
%
%   Background: Based on p34.m from [1].
%
%   Reference: [1] L. N. Trefethen, Spectral methods in MATLAB. Society
%              for Industrial and Applied Mathematics, 2000.
%              doi: 10.1137/1.9780898719598.
%
%   Part of the PPR repository.
%%
if nargin < 4
    if nargin < 3
        if nargin < 2
            if nargin < 1
                n = 14;
            end
            degree = 6;
        end
        eps = 0.5;
    end
    r = 14;
end

% clear; close all; degree = 6; n = 14; r = 14;
plots = 'animation';

fprintf('Running Example 9 with ε=%.1f\n',eps)
%% Get dynamics and define control problem
[f, B, Q, ~, y] = getSystem9Neumann(eps, n);
R = 1e-1; m = size(B,2);

%% Construct controllers
% Open-loop (uncontrolled) controller
uUnc = @(z) zeros(m,1);

% PPR Controller
fprintf("Computing ppr() solution, n=%i, r=%i, d=%i ... ",n,r,degree); tic
options = struct; options.verbose = false; options.r = r; options.h = B.';
[~, GainsPPR, options] = ppr(f, B, Q, R, degree, options);
fprintf("completed in %2.2f seconds. \n", toc)

uPPR = @(z) (kronPolyEval(GainsPPR, z));

% LQR Controller (first term in PPR controller)
uLQR = @(z) (kronPolyEval(GainsPPR, z, 1));

% SDRE Controller
uSDRE = @(z) sdre(@(y)(f{1}+diag(y.^2)),@(y)(B),Q,R,z);

% TT-HJB Controller
addpath(genpath('../TT-Toolbox/'),genpath('../tamen/'),genpath('../TT-HJB/'))
rng("default")
umax = inf;nv = 5;av = 3;tol = 1e-3;mu = 100;
% Prepare handles for ODE component functions for HJB
ffun = @(x,i)x*f{1}(i,:).' - x(:,i).^3; % nonlinear Dynamics
gfun = @(x,i)B(i)*ones(size(x,1),1); % Actuator (here just const vector)
lfun = @(x)sum((x.^2).*diag(Q).', 2);
% Call Dolgov code
fprintf("Computing TT-HJB solution ... "); figure; tic
[printout, V] = evalc('hjb_leg(size(f{1},1), nv, av, ffun, gfun, lfun, R, tol, mu, [], umax)'); % Call with evalc to keep command window clean
uTTHJB = @(z)controlfun_leg(z,-av,av,core2cell(V),gfun,R,umax);
fprintf("completed in %2.2f seconds. \n", toc)
% close % close convergence plot
% uTTHJB =  @(z) uLQR(z.'); % For debugging my stuff without running tthjb for so long

%% Simulate closed-loop systems
% Construct original system dynamics
FofXU = @(v,u) (f{1}*v - v.^3 + B*u); % -v.^3 is equivalent to f{3}*x^(3), just faster for simulation purposes

switch plots
    case 'animation'
        t = [0,1e7];
    case 'figs'
        tmax = 5; dt = 0.001; t = 0:dt:tmax; % Specify time vector to accurately approximate cost function integral
end

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-20);

vOffsets = [0, 0.25, 0.5, 1, 1.5, 2];
for j=1:6
    clear UxPPR UxTTHJB UxSDRE
    v0 = vOffsets(j) + cos(2*pi*y).*cos(pi*y); % Modified initial condition

    % Simulate using ode solver (more efficient than forward euler)
    [tUnc, Xunc] = ode15s(@(t, v) FofXU(v,uUnc(v)),t, v0, opts);
    [tLQR, XLQR] = ode15s(@(t, v) FofXU(v,uLQR(v)),t, v0, opts);
    [tPPR, XPPR] = ode15s(@(t, v) FofXU(v,uPPR(v)),t, v0, opts);
    [tSDRE, XSDRE] = ode15s(@(t, v) FofXU(v,uSDRE(v)),t, v0, opts);
    [tTTHJB, XTTHJB] = ode15s(@(t, v) FofXU(v,uTTHJB(v.')),t, v0, opts);


    % Compute performance index (cost)
    UxUnc = uUnc(Xunc.').'; UxLQR = uLQR(XLQR.').';
    for i=1:length(tPPR); UxPPR(i,1) = uPPR(XPPR(i,:).'); end
    for i=1:length(tSDRE); UxSDRE(i,1) = uSDRE(XSDRE(i,:).'); end
    for i=1:length(tTTHJB); UxTTHJB(i,1) = uTTHJB(XTTHJB(i,:)); end
    costUnc(j) = trapz(tUnc, sum((Xunc.^2).*diag(Q).', 2));
    costLQR(j) = trapz(tLQR, sum((XLQR.^2).*diag(Q).', 2) + R*UxLQR.^2);
    costPPR(j) = trapz(tPPR, sum((XPPR.^2).*diag(Q).', 2) + R*UxPPR.^2);
    costSDRE(j) = trapz(tSDRE, sum((XSDRE.^2).*diag(Q).', 2) + R*UxSDRE.^2);
    costTTHJB(j) = trapz(tTTHJB, sum((XTTHJB.^2).*diag(Q).', 2) + R*UxTTHJB.^2);

    % fprintf("  Controller costs:")
    % fprintf("\n    Uncontrolled    & %3.3f         ",costUnc(j))
    % fprintf("\n    LQR             & %3.3f         ",costLQR(j))
    % fprintf("\n    PPR             & %3.3f         ",costPPR(j))
    % fprintf("\n    SDRE            & %3.3f         ",costSDRE(j))
    % fprintf("\n    TT-HJB          & %3.3f         ",costTTHJB(j))

    switch plots
        case 'animation'
            figure
            title('Closed-loop animation')
            for j = 1:length(tUnc)
                tt = tUnc(j);
                hold off
                [~,i] = min(abs(tUnc-tt));
                plot(Xunc(i,:));
                hold on
                [~,i] = min(abs(tLQR-tt));
                plot(XLQR(i,:));
                [~,i] = min(abs(tPPR-tt));
                plot(XPPR(i,:));
                [~,i] = min(abs(tSDRE-tt));
                plot(XSDRE(i,:)); ylim([-1 1]);
                [~,i] = min(abs(tTTHJB-tt));
                plot(XTTHJB(i,:)); ylim([-1 1]);
                ylim([-3 3]); xlim([0 14])
                legend('Uncontrolled','LQR','PPR','SDRE','TT-HJB')
                pause(0.01)
            end
        case 'figs'
            % Plots just for checking results, not for the paper
            plotIndices = round(linspace(0,1,51).^3*5000+1);
            % plotIndices = linspace(1,5001,51);
            plotT = t(1:100:end); nplots = length(plotT);
            xx = -1:.025:1; plotdata = [zeros(nplots,length(xx))];
            X = Xunc(plotIndices,:);
            for i=1:length(plotT)
                plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
            end
            figure, subplot('position',[.1 .4 .8 .5])
            mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -2.05 3.05]),
            view(-60,55), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
            title("Open-loop"); drawnow

            X = XLQR(plotIndices,:);
            for i=1:length(plotT)
                plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
            end
            figure, subplot('position',[.1 .4 .8 .5])
            mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -2.05 3.05]),
            view(-60,55), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
            title("LQR"); drawnow

            X = XPPR(plotIndices,:);
            for i=1:length(plotT)
                plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
            end
            figure, subplot('position',[.1 .4 .8 .5])
            mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -2.05 3.05]),
            view(-60,55), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
            title("PPR"); drawnow

            X = XTTHJB(plotIndices,:);
            for i=1:length(plotT)
                plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
            end
            figure, subplot('position',[.1 .4 .8 .5])
            mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -2.05 3.05]),
            view(-60,55), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
            title("TT-HJB"); drawnow
    end
end

fprintf('\n# Table II Data (Allen-Cahn, Neumann BCs)\n');
fprintf('# Control costs for different initial condition offsets\n');
fprintf("      Controller    &     v0=%2.2f      &     v0=%2.2f      &     v0=%2.2f      &     v0=%2.2f      &     v0=%2.2f      &     v0=%2.2f     \n",vOffsets)
    fprintf("     %s   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n", 'Uncontrolled', costUnc)
    fprintf("     %s   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n", '     LQR    ', costLQR)
    fprintf("     %s   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n", '     PPR    ', costPPR)
    fprintf("     %s   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n", '     SDRE   ', costSDRE)
    fprintf("     %s   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n\n", '     TTHJB  ', costTTHJB)
end
