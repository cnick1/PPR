function runExample10()
%runExample10 Runs the 1D transcritical bifurcation example.
%
%   Usage:  runExample10()
%
%   Background:
%
%   Reference:
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


%% See if closed-loop simulations match expectation
return; % Comment this out to run the rest of the script
uSDRE = @(z) sdre(@(y)(f{1}+f{2}*y+f{3}*y^2),@(y)(g{1}),Q2+diag(z.^2),R,z);
controllers = {uOpenLoop, uLQR, uSDRE, uPPR};
controllerNames = {'Uncontrolled', 'LQR         ', 'SDRE        ', 'PPR         '};
figure('Position',[113 277.6667 1.5807e+03 482.6667]);
for idx = 1:length(controllers)
    u = controllers{idx};

    % Simulate using ode solver (more efficient than forward euler)
    subplot(2,2,idx); title(sprintf("Controller %s ",controllerNames{idx})); hold on;
    v0s = -3:0.1:3;
    for j=1:length(v0s)
        v0 = v0s(j);
        try
            [t, X] = ode23s(@(t, x) FofXU(x,u(x)),[0,5], v0);

            % Compute performance index (cost)
            Lagrangian = zeros(size(t));
            for i=1:length(t)
                xbar = X(i,:).'; Ux = u(xbar);
                Lagrangian(i) = 1/2*(Q2*xbar^2 + R*Ux^2 + q{3}*xbar^3); % hardcoded v.^4 instead of Q4 for speed
            end
            performanceIndex(idx,j) = trapz(t, Lagrangian);

            plot(X,t)
            xlim([-3, 3]); drawnow
        catch
            performanceIndex(idx,j) = inf;
        end
    end
end

figure('Position',[113 277.6667 1.5807e+03 482.6667]);
for idx = 1:length(controllers)
    u = controllers{idx};
    % Simulate using ode solver (more efficient than forward euler)
    subplot(2,2,idx); title(sprintf("Controller %s ",controllerNames{idx})); hold on;
    v0s = -3:0.01:3;
    for j=1:length(v0s)
        v0 = v0s(j);
        plot(v0,u(v0),'k.')
        xlim([-3, 3]); drawnow
    end
end

end
