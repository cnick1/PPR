function runExample9_differentICs(n, degree, r)
%runExample9_differentICs Runs the Allen-Cahn example.
%
%   Usage:  runExample9(n,degree,r)
%
%   Inputs: n      - desired state dimension
%           degree - desired polynomial degree of value function to compute
%           r      - ROM dimension; if r=n, no MOR is performed
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
            n = 33;
        end
        degree = 4;
    end
    r = 10;
end

fprintf('Running Example 9, Allen-Cahn example with Dirichlet BCs, for different initial conditions \n')

%% Construct controller
y0 = .5; % Desired interface location
eps = 0.01;

% Get system expanded about vref, reference configuration (@ origin) -> v = v+vref
[f, B, ~, D, y, vref] = getSystem9(eps, n-1, y0);

B = B(:,linspace(1,n,5)); B(:,[1 5]) = []; m = size(B,2);
Q2 = 0.1; Q3 = sparse(n^3,1) ; Q4 = sparse(linspace(1,n^4,n),1,1);
q = {[],Q2,Q3,Q4}; R = 1;

% Full PPR solution (LQR is just the first term)
fprintf("Computing ppr() solution, n=%i, d=%i ... ",n,degree); tic
options = struct; options.verbose = false; options.h = B.';
[~, GainsPPR, options] = ppr(f, B, q, R, degree, options);
fprintf("completed in %2.2f seconds. \n", toc)

% Reduced PPR Solution
fprintf("Computing ppr() solution, n=%i, r=%i, d=%i ... ",n,r,degree); tic
options = struct; options.verbose = false; options.r = r; options.h = B.';
[~, GainsPPR_reduced, options] = ppr(f, B, q, R, degree, options);
fprintf("completed in %2.2f seconds. \n", toc)

% LPR Solution
[~, GainsLPR] = ppr(f(1), B, q, R, degree, options);

% Construct control laws
uOpenLoop = @(z) zeros(m,1);
uLQR = @(z) (kronPolyEval(GainsPPR, z, 1));
uLPR= @(z) (kronPolyEval(GainsLPR, z));
uSDRE = @(z) sdre(@(y)(f{1}-3*diag(vref).*diag(y)-diag(y.^2)),@(y)(B),Q2+diag(z.^2),R,z);
uPPR = @(z) (kronPolyEval(GainsPPR, z));
uPPR_reduced = @(z) (kronPolyEval(GainsPPR_reduced, z));

controllers = {uOpenLoop, uLQR, uLPR, uSDRE, uPPR, uPPR_reduced};
controllerNames = {'Uncontrolled', 'LQR         ', 'LPR         ', 'SDRE        ', 'PPR         ', 'PPR reduced '};

%% Simulate closed-loop systems
% Construct original system dynamics
D2 = D^2; D2([1 n],:) = zeros(2,n); % For boundary conditions
FofXU = @(v,u) (eps*D2*v + v - v.^3 + B*u);

tmax = 100; dt = .02; t = 0:dt:tmax; % Specify time vector for plotting

x0s = -.75:.25:.75; performanceIndex=zeros(6,length(x0s));
for j = 1:length(x0s)
    x0 = x0s(j);

    % Initial condition with x0 interface location
    v0 = tanh((y-(x0))/sqrt(2*eps));

    for idx = 1:length(controllers)
        u = controllers{idx};

        % Simulate using ode solver (more efficient than forward euler)
        try
            [t, X] = ode23s(@(t, v) FofXU(v,u(v-vref)),t, v0);

            % Compute performance index (cost)
            Lagrangian = zeros(size(t));
            for i=1:length(t)
                xbar = X(i,:).' - vref; Ux = u(xbar);
                Lagrangian(i) = 1/2*(xbar.'*Q2*xbar + Ux.'*R*Ux + sum(xbar.^4)); % hardcoded v.^4 instead of Q4 for speed
            end
            performanceIndex(idx,j) = trapz(t, Lagrangian);

            % Plots just for checking results, not for the paper
            plotT = t(1:100:end); X = X(1:100:end,:); nplots = length(plotT);
            xx = -1:.025:1; plotdata = [zeros(nplots,length(xx))];
            for i=1:length(plotT)
                plotdata(i,:) = polyval(polyfit(y,X(i,:),20),xx);
            end
            figure, subplot('position',[.1 .4 .8 .5])
            mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -1.05 1.05]),
            view(-60,55), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
            title(sprintf("Controller %s  (x0=%2.2f)",controllerNames{idx},x0)); drawnow
        catch
            performanceIndex(idx,j) = inf;
        end
    end
end

fprintf('\n# Table II Data (Allen-Cahn, Dirichlet BCs)\n');
fprintf('# Control costs for different initial conditions\n');
fprintf("      Controller    &   x0=%5.2f   &   x0=%5.2f   &   x0=%5.2f   &   x0=%5.2f   &   x0=%5.2f   &   x0=%5.2f   &   x0=%5.2f  \n",x0s)
for idx = 1:length(controllers)
    fprintf("     %s   &  %9.3f   &  %9.3f   &  %9.3f   &  %9.3f   &  %9.3f   &  %9.3f   &  %9.3f       \n",controllerNames{idx},performanceIndex(idx,:))
end

% exportgraphics(figure(6),sprintf('plots/example9_mor_n%i_r%i_d%i.pdf',n,r,degree), 'ContentType', 'vector')
% exportgraphics(figure(2),sprintf('plots/example9_mor_n%i_d%i.pdf',n,degree), 'ContentType', 'vector')
% exportgraphics(figure(1),sprintf('plots/example9_mor_openloop_n%i.pdf',n), 'ContentType', 'vector')

end
