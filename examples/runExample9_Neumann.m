function runExample9_Neumann(n, degree, eps, vOffset, r)
%runExample9_Neumann Runs the Allen-Cahn example with Nuemman BCs.
%
%   Usage:  runExample9_Neumann(n,degree,eps,vOffset,r)
%
%   Inputs: n       - desired state dimension
%           degree  - desired polynomial degree of value function to compute
%           eps     - diffusion coefficient
%           vOffset - shift of initial condition from the origin
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
if nargin < 5
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
        vOffset = 0.1;
    end
    r = 14;
end

% clear; close all; degree = 6; n = 14; r = 14;
plots = 'animation';

fprintf('Running Example 9 with ε=%.1f, v0=%.1f\n',eps,vOffset)

%% Construct controller
% Get system expanded about vref, reference configuration (@ origin) -> v = v+vref
[f, B, Q, ~, y] = getSystem9Neumann(eps, n);
R = 1e-1; m = size(B,2);

% Compute PPR solution (LQR is just the first term)
fprintf("Computing ppr() solution, n=%i, r=%i, d=%i ... ",n,r,degree); tic
options = struct; options.verbose = false; options.r = r; options.h = B.';
[~, GainsPPR, options] = ppr(f, B, Q, R, degree, options);
fprintf("completed in %2.2f seconds. \n", toc)

uUnc = @(z) zeros(m,1);
uLQR = @(z) (kronPolyEval(GainsPPR, z, 1));
uPPR = @(z) (kronPolyEval(GainsPPR, z));
uSDRE = @(z) sdre(@(y)(f{1}+diag(y.^2)),@(y)(B),Q,R,z);

%% Solve with TT-HJB
addpath(genpath('../TT-Toolbox/'),genpath('../tamen/'),genpath('../TT-HJB/'))
rng("default")
umax = inf;nv = 5;av = 3;tol = 1e-3;mu = 100;
% Prepare handles for ODE component functions for HJB
ffun = @(x,i)x*f{1}(i,:).' - x(:,i).^3; % nonlinear Dynamics
gfun = @(x,i)B(i)*ones(size(x,1),1); % Actuator (here just const vector)
lfun = @(x)sum((x.^2).*diag(Q).', 2);

% Solve HJB-controlled ODE
fprintf("Computing TT-HJB solution ... "); figure; tic
% [printout, V] = evalc('hjb_leg(size(f{1},1), nv, av, ffun, gfun, lfun, R, tol, mu, [], umax)');
% uTTHJB = @(z)controlfun_leg(z,-av,av,core2cell(V),gfun,R,umax);
fprintf("completed in %2.2f seconds. \n", toc)
close
uTTHJB =  @(z) uLQR(z.'); % For debugging my stuff without running tthjb for so long
%
%% Simulate closed-loop systems
% Construct original system dynamics
FofXU = @(v,u) (f{1}*v - v.^3 + B*u); % -v.^3 is equivalent to f{3}*x^(3), just faster for simulation purposes

v0 = vOffset + cos(2*pi*y).*cos(pi*y); % Modified initial condition
% v0 = .53*y + .47*sin(-1.5*pi*y); % Initial condition from Trefethen

switch plots
    case 'animation'
        t = [0,1e6];
    case 'figs'
        tmax = 5; dt = 0.001; t = 0:dt:tmax; % Specify time vector to accurately approximate cost function integral
end

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-20);

% Simulate using ode solver (more efficient than forward euler)
clear UxPPR UxTTHJB
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
LagrangianUnc = sum((Xunc.^2).*diag(Q).', 2);
LagrangianLQR = sum((XLQR.^2).*diag(Q).', 2) + R*UxLQR.^2;
LagrangianPPR = sum((XPPR.^2).*diag(Q).', 2) + R*UxPPR.^2;
LagrangianSDRE = sum((XSDRE.^2).*diag(Q).', 2) + R*UxSDRE.^2;
LagrangianTTHJB = sum((XTTHJB.^2).*diag(Q).', 2) + R*UxTTHJB.^2;

fprintf("  Controller costs:")
fprintf("\n    Uncontrolled    & %3.3f         ",trapz(tUnc, LagrangianUnc))
fprintf("\n    LQR             & %3.3f         ",trapz(tLQR, LagrangianLQR))
fprintf("\n    PPR             & %3.3f         ",trapz(tPPR, LagrangianPPR))
fprintf("\n    SDRE            & %3.3f         ",trapz(tSDRE, LagrangianSDRE))
fprintf("\n    TT-HJB          & %3.3f         ",trapz(tTTHJB, LagrangianTTHJB))

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

fprintf('\n\n')

end
