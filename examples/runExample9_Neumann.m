function runExample9_Neumann(n, degree, eps)
%runExample9_Neumann Runs the Allen-Cahn example with Nuemman BCs.
%
%   Usage:  runExample9_Neumann(n,degree,eps)
%
%   Inputs: n       - desired state dimension
%           degree  - desired polynomial degree of value function to compute
%           eps     - diffusion coefficient
%
%   Background: Based on p34.m from [1].
%
%   Reference: [1] L. N. Trefethen, Spectral methods in MATLAB. Society
%              for Industrial and Applied Mathematics, 2000.
%              doi: 10.1137/1.9780898719598.
%
%   Part of the PPR repository.
%%
if nargin < 3
    if nargin < 2
        if nargin < 1
            n = 14;
        end
        degree = 6;
    end
    eps = 0.5;
end
% Create a vec function for readability
vec = @(X) X(:);

% clear; close all; degree = 6; n = 14;
r = 3;
plots = 'figs';

fprintf('\nRunning Example 9, Allen-Cahn example with Nuemman BCs, for ε=%.1f\n',eps)

%% Get dynamics and define control problem
[f, B, Q, ~, y] = getSystem9Neumann(eps, n);
R = 1e-1; m = size(B,2);

%% Construct controllers
% Open-loop (uncontrolled) controller
uUnc = @(z) zeros(m,1);

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
options = struct; options.verbose = false; options.r = r; options.h = B.';
[~, GainsPPR_reduced, options] = ppr(f, B, Q, R, 6, options);
fprintf("completed in %2.2f seconds. \n", toc)

uPPR_reduced = @(z) (kronPolyEval(GainsPPR_reduced, z));


% Tuned Reduced PPR Controller
fprintf(" Computing tuned reduced ppr() solution, n=%i, r=%i, d=%i ... ",n,r,4); tic
options = struct; options.verbose = false; options.r = r; options.h = B.';
[~, GainsPPR_tuned_reduced, options] = ppr(f, B, {0,vec(Q),0,0.25}, {R,0.009,0}, 4, options);
fprintf("completed in %2.2f seconds. \n", toc)

uPPR_tuned_reduced = @(z) (kronPolyEval(GainsPPR_tuned_reduced, z));

% LQR Controller (first term in PPR controller)
uLQR = @(z) (kronPolyEval(GainsPPR, z, 1));

% SDRE Controller
uSDRE = @(z) sdre(@(y)(f{1} - diag(y.^2)),@(y)(B),Q,R,z);
% uSDRE =  @(z) uLQR(z); % For debugging my stuff without running sdre for so long

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

switch plots
    case 'animation'
        t = [0,1e7];
    case 'figs'
        tmax = 5; dt = 0.001; t = 0:dt:tmax; % Specify time vector for plotting
end

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-20);

wOffsets = [0.5, 1, 1.5, 2];
% vOffsets = [0, 0.25, 0.5, 1, 1.5, 2];
for j=1:length(wOffsets)
    clear UxPPR UxPPR_tuned UxPPR_tuned_reduced UxPPR_reduced UxTTHJB UxSDRE
    w0 = wOffsets(j) + cos(2*pi*y).*cos(pi*y); % Modified initial condition
    
    % Simulate using ode solver (more efficient than forward euler) and
    % compute performance indexes
    [tUnc, Xunc] = ode15s(@(t, v) FofXU(v,uUnc(v)),t, w0, opts);
    costUnc(j) = trapz(tUnc, sum((Xunc.^2).*diag(Q).', 2));
    
    [tLQR, XLQR] = ode15s(@(t, v) FofXU(v,uLQR(v)),t, w0, opts);
    UxLQR = uLQR(XLQR.').';
    costLQR(j) = trapz(tLQR, sum((XLQR.^2).*diag(Q).', 2) + R*UxLQR.^2);
    
    [tPPR, XPPR] = ode15s(@(t, v) FofXU(v,uPPR(v)),t, w0, opts);
    for i=1:length(tPPR); UxPPR(i,1) = uPPR(XPPR(i,:).'); end
    costPPR(j) = trapz(tPPR, sum((XPPR.^2).*diag(Q).', 2) + R*UxPPR.^2);
    
    [tPPR_tuned, XPPR_tuned] = ode15s(@(t, v) FofXU(v,uPPR_tuned(v)),t, w0, opts);
    for i=1:length(tPPR_tuned); UxPPR_tuned(i,1) = uPPR_tuned(XPPR_tuned(i,:).'); end
    costPPR_tuned(j) = trapz(tPPR_tuned, sum((XPPR_tuned.^2).*diag(Q).', 2) + R*UxPPR_tuned.^2);
    
    [tPPR_tuned_reduced, XPPR_tuned_reduced] = ode15s(@(t, v) FofXU(v,uPPR_tuned_reduced(v)),t, w0, opts);
    for i=1:length(tPPR_tuned_reduced); UxPPR_tuned_reduced(i,1) = uPPR_tuned_reduced(XPPR_tuned_reduced(i,:).'); end
    costPPR_tuned_reduced(j) = trapz(tPPR_tuned_reduced, sum((XPPR_tuned_reduced.^2).*diag(Q).', 2) + R*UxPPR_tuned_reduced.^2);
    
    [tPPR_reduced, XPPR_reduced] = ode15s(@(t, v) FofXU(v,uPPR_reduced(v)),t, w0, opts);
    for i=1:length(tPPR_reduced); UxPPR_reduced(i,1) = uPPR_reduced(XPPR_reduced(i,:).'); end
    costPPR_reduced(j) = trapz(tPPR_reduced, sum((XPPR_reduced.^2).*diag(Q).', 2) + R*UxPPR_reduced.^2);
    
    try
        [tSDRE, XSDRE] = ode15s(@(t, v) FofXU(v,uSDRE(v)),t, w0, opts);
        for i=1:length(tSDRE); UxSDRE(i,1) = uSDRE(XSDRE(i,:).'); end
        costSDRE(j) = trapz(tSDRE, sum((XSDRE.^2).*diag(Q).', 2) + R*UxSDRE.^2);
    catch
        costSDRE(j) = inf;
    end
    
    [tTTHJB, XTTHJB] = ode15s(@(t, v) FofXU(v,uTTHJB(v.')),t, w0, opts);
    for i=1:length(tTTHJB); UxTTHJB(i,1) = uTTHJB(XTTHJB(i,:)); end
    costTTHJB(j) = trapz(tTTHJB, sum((XTTHJB.^2).*diag(Q).', 2) + R*UxTTHJB.^2);
    
    switch plots
        case 'animation'
            fig2 = figure('Position',[664.3333 1.5323e+03 1.0133e+03 557.3333]);
            plot(tUnc,tUnc*0)
            hold on;
            plot(tLQR,UxLQR)
            plot(tPPR,UxPPR)
            plot(tPPR_tuned,UxPPR_tuned)
            plot(tPPR_reduced,UxPPR_reduced)
            plot(tPPR_tuned_reduced,UxPPR_tuned_reduced)
            plot(tSDRE,UxSDRE)
            plot(tTTHJB,UxTTHJB)
            legend('Uncontrolled','LQR','PPR','PPR tuned','PPR reduced','PPR tuned reduced','SDRE','TTHJB')
            xlim([0 5])
            drawnow
            
            figure
            title('Closed-loop animation')
            for j2 = 1:length(tUnc)
                tt = tUnc(j2);
                hold off
                [~,i] = min(abs(tUnc-tt));
                plot(Xunc(i,:));
                hold on
                [~,i] = min(abs(tLQR-tt));
                plot(XLQR(i,:));
                [~,i] = min(abs(tPPR-tt));
                plot(XPPR(i,:));
                [~,i] = min(abs(tPPR_tuned-tt));
                plot(XPPR_tuned(i,:));
                [~,i] = min(abs(tPPR_reduced-tt));
                plot(XPPR_reduced(i,:));
                [~,i] = min(abs(tPPR_tuned_reduced-tt));
                plot(XPPR_tuned_reduced(i,:));
                [~,i] = min(abs(tSDRE-tt));
                plot(XSDRE(i,:)); ylim([-1 1]);
                [~,i] = min(abs(tTTHJB-tt));
                plot(XTTHJB(i,:)); ylim([-1 1]);
                ylim([-3 3]); xlim([0 14])
                legend('Uncontrolled','LQR','PPR','PPR tuned','PPR reduced','PPR tuned reduced','SDRE','TTHJB')
                pause(0.005)
            end
        case 'figs'
            fig2 = figure;
            plot(tUnc,tUnc*0)
            hold on;
            plot(tLQR,UxLQR)
            plot(tPPR,UxPPR)
            plot(tPPR_reduced,UxPPR_reduced)
            plot(tPPR_tuned_reduced,UxPPR_tuned_reduced)
            % plot(tSDRE,UxSDRE)
            % plot(tTTHJB,UxTTHJB)
            legend('Uncontrolled','LQR','PPR','PPR reduced','PPR tuned reduced')%,'SDRE','TT-HJB')
            %%
            figure('Position',[97 164.3333 1558 577.3333])
            % Plots just for checking results, not for the paper
            plotIndices = round(linspace(0,1,51).^3*5000+1);
            % plotIndices = linspace(1,5001,51);
            plotT = t(1:100:end); nplots = length(plotT);
            xx = -1:.025:1; plotdata = [zeros(nplots,length(xx))];
            X = Xunc(plotIndices,:);
            for i=1:length(plotT)
                plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
            end
            subplot(2,4,1)
            mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -0.5 3.05]),
            view(145,35), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
            title("Open-loop"); drawnow
            
            X = XLQR(plotIndices,:);
            for i=1:length(plotT)
                plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
            end
            subplot(2,4,2)
            mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -0.5 3.05]),
            view(145,35), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
            title("LQR"); drawnow
            
            X = XSDRE(plotIndices,:);
            for i=1:length(plotT)
                plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
            end
            subplot(2,4,3)
            mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -0.5 3.05]),
            view(145,35), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
            title("SDRE"); drawnow
            
            X = XTTHJB(plotIndices,:);
            for i=1:length(plotT)
                plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
            end
            subplot(2,4,4)
            mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -0.5 3.05]),
            view(145,35), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
            title("TT-HJB"); drawnow
            
            
            X = XPPR(plotIndices,:);
            for i=1:length(plotT)
                plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
            end
            subplot(2,4,5)
            mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -0.5 3.05]),
            view(145,35), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
            title("PPR"); drawnow
            
            X = XPPR_reduced(plotIndices,:);
            for i=1:length(plotT)
                plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
            end
            subplot(2,4,6)
            mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -0.5 3.05]),
            view(145,35), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
            title("PPR reduced"); drawnow
            
            X = XPPR_tuned(plotIndices,:);
            for i=1:length(plotT)
                plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
            end
            subplot(2,4,7)
            mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -0.5 3.05]),
            view(145,35), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
            title("PPR tuned"); drawnow
            
            X = XPPR_tuned_reduced(plotIndices,:);
            for i=1:length(plotT)
                plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
            end
            subplot(2,4,8)
            mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -0.5 3.05]),
            view(145,35), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
            title("PPR tuned reduced"); drawnow
            
            exportgraphics(gcf,sprintf('plots/example9_neumann_v0%2.1f.pdf',wOffsets(j)), 'ContentType', 'vector')
            
            %%
    end
end


fprintf('\n# Table III Data (Allen-Cahn, Neumann BCs)\n');
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


% fprintf('\n# Table III Data (Allen-Cahn, Neumann BCs)\n');
% fprintf('# Control costs for different initial condition offsets\n');
% fprintf("      Controller    &     v0=%2.2f      &     v0=%2.2f      &     v0=%2.2f      &     v0=%2.2f      &     v0=%2.2f      &     v0=%2.2f     \n",wOffsets)
% fprintf("     %s   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n", 'Uncontrolled', costUnc)
% fprintf("     %s   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n", '     LQR    ', costLQR)
% fprintf("     %s   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n", '     SDRE   ', costSDRE)
% fprintf("     %s   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n", '     PPR    ', costPPR)
% fprintf("     %s   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n", ' PPR reduced', costPPR_reduced)
% fprintf("     %s   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n", ' PPR tuned  ', costPPR_tuned)
% fprintf(" %s  &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n", 'PPR tuned reduced', costPPR_tuned_reduced)
% fprintf("     %s   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f   &  %13.3f       \n\n", '     TTHJB  ', costTTHJB)
end


% for i=1:8
%     subplot(2,4,i)
%     axis([-1 1 0 tmax -1 3.05])
% end

