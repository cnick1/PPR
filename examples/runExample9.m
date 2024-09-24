function runExample9(n, degree)
%runExample9 Runs the Allen-Cahn example.
%
%   Usage:  runExample9
%
%   Inputs: n      - desired state dimension
%           degree - desired polynomial degree of value function to compute
%   Outputs:
%
%   Background: Based on p34.m from [1].
%
%   Reference: [1] L. N. Trefethen, Spectral methods in MATLAB. Society
%              for Industrial and Applied Mathematics, 2000.
%              doi: 10.1137/1.9780898719598.
%
%   Part of the NLbalancing repository.
%%

if nargin < 2
    degree = 4;
    % degree = 2;
    if nargin < 1
        % n = 129;
        n = 33;
    end
end

fprintf('Running Example 9\n')

%% Construct controller
y0 = .5; % Desired interface location

for eps = [0.01 0.0075 0.005]
    % Get system expanded about vref, reference configuration (@ origin) -> v = v+vref
    [f, B, ~, D, y, vref] = getSystem9(eps, n-1, y0);

    B = B(:,linspace(1,n,5)); B(:,[1 5]) = []; m = size(B,2);
    Q2 = 0.1; Q3 = sparse(n^3,1) ; Q4 = sparse(linspace(1,n^4,n),1,4);
    q = {[],Q2,Q3,Q4}; R = 1;

    % Compute PPR solution (LQR is just the first term)
    fprintf("Computing ppr() solution, n=%i, d=%i ... ",n,degree)
    options.verbose = false; options.skipGains = false;
    [~, GainsPPR] = ppr(f, B, q, R, degree, options);
    fprintf("completed.\n")

    uOpenLoop = @(z) zeros(m,1);
    uLQR = @(z) (kronPolyEval(GainsPPR, z, 1));
    uPPR_quad = @(z) (kronPolyEval(GainsPPR, z, 2));
    uPPR = @(z) (kronPolyEval(GainsPPR, z));
    controllers = {uOpenLoop, uLQR, uPPR_quad, uPPR};

    %% Simulate closed-loop systems
    % Construct original system dynamics
    D2 = D^2; D2([1 n],:) = zeros(2,n); % For boundary conditions
    FofXU = @(v,u) (eps*D2*v + v - v.^3 + B*u);

    % Initial condition from Trefethen
    v0 = .53*y + .47*sin(-1.5*pi*y);

    tmax = 1000; dt = .2; t = 0:dt:tmax; % Specify time vector to accurately approximate cost function integral

    fprintf("  Controller & Cost (eps=%2.4f)    ",eps)
    for idx = 1:4
        u = controllers{idx};

        % Simulate using ode solver (more efficient than forward euler)
        [t, X] = ode23s(@(t, v) FofXU(v,u(v-vref)),t, v0);

        % Compute performance index (cost)
        Lagrangian = zeros(size(t));
        for i=1:length(t)
            xbar = X(i,:).' - vref; Ux = u(xbar);
            Lagrangian(i) = 1/2*(xbar.'*Q2*xbar + Ux.'*R*Ux + 4*sum(xbar.^4)); % hardcoded v.^4 instead of Q4 for speed
        end
        performanceIndex = trapz(t, Lagrangian);
        fprintf("\n        %i    & %3.3f         ",idx-1,performanceIndex)

        % Plots just for checking results, not for the paper
        plotT = t(1:100:end); X = X(1:100:end,:); nplots = length(plotT);
        xx = -1:.025:1; plotdata = [zeros(nplots,length(xx))];
        for i=1:length(plotT)
            plotdata(i,:) = polyval(polyfit(y,X(i,:),20),xx);
        end
        figure, subplot('position',[.1 .4 .8 .5])
        mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -1.05 1.05]),
        view(-60,55), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
        title(sprintf("Controller %i  (eps=%2.4f)",idx-1,eps)); drawnow
    end
    fprintf('\n')
end

end
