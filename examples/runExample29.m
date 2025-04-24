function runExample29(numElements, r, descriptor)
%runExample29 Runs the 2D unsteady nonlinear heat equation FEM control example.
%
%   Usage:  runExample29(numElements)
%
%   Inputs:
%       numElements - number of finite elements in each direction
%       r           - reduced-order dimension
%       descriptor  - boolean, option to leave dynamics in generalized form
%                     (default=true)
%
%   Description: The model, inspired by [1], describes an unstable heat
%   equation problem modeling heat generated in a resistive electrical
%   material. The resulting model takes the form
%
%     uₜ(x,y,t) = uₓₓ(x,y,t) + uᵧᵧ(x,y,t) + λ u(x,y,t)
%     uᵧ(x,0,t) = u₁(t)  (Neumann control input on side AB)
%     uₓ(1,y,t) = u₂(t)  (Neumann control input on side BC)
%     uᵧ(x,1,t) = u₃(t)  (Neumann control input on side CD)
%     uₓ(0,y,t) = u₄(t)  (Neumann control input on side DA)
%
%   where one end is insulated and one end is subject to Neumann boundary
%   control. The FEM model can be written as
%
%     M ẋ + K₁ x + K₃ (x⊗x⊗x) = B u,
%     y = C x.
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
fprintf('Running Example 29\n')
if nargin < 3
    descriptor = true;
    if nargin < 2
        if nargin < 1
            numElements = 8;
        end
        r = (numElements+1)^2;
    end
end

% Get dynamics
nx = numElements+1; ny = nx;
n = nx*ny; m = 1;
fprintf(" Forming FEM model, n=%i ... ",n); tic
[E, f, g, ~, xyg] = getSystem29(numElements,.75,1,-1);
fprintf("completed in %2.2f seconds. \n", toc)
% Insulate all sides, use control only on side CD
g{1} = g{1}(:,3);
G = @(x) g{1};

% Pass sparse operators for ode simulation purposes
Mchol = chol(E).'; % Use Cholesky factor for inverting
FofXU = @(x,u) Mchol.'\(Mchol\(kronPolyEval({sparse(f{1}),f{2},f{3}},x) + g{1} * u));

if descriptor
    % Case 1: Descriptor form (sparse)
    options.E = E;
    fprintf(" Computing ppr() solution in descriptor form, n=%i, r=%i, d=%i ... ",n,r,4); tic
else
    % Case 2: Standard form (dense)
    f{1} = Mchol.'\(Mchol\f{1});
    f{3} = Mchol.'\(Mchol\f{3});
    g{1} = Mchol.'\(Mchol\g{1});
    options.E = []; % Case 3: Standard form solved with E=I
    fprintf(" Computing ppr() solution in standard form, n=%i, r=%i, d=%i ... ",n,r,4); tic
end

% Get value function/controller
q = .1; R = 1; degree = 4; options.verbose = true; options.r = r;
[~, K] = ppr(f, g, q, R, degree, options);
fprintf("completed in %2.2f seconds. \n", toc)

uLQR = @(x) kronPolyEval(K, x, 1);
uPPR = @(x) kronPolyEval(K, x, degree-1);

%% Simulate closed-loop system
X = reshape(xyg(:,1),nx,ny);
Y = reshape(xyg(:,2),nx,ny);
x0 = .25*(sin(4*pi*X) + cos(3*pi*Y)) + .1;
x0 = x0(:);
tmax = 5; t = 0:0.2:tmax; % specify for plotting
opts=odeset('OutputFcn',@odeprog);

fprintf(" - Simulating open-loop dynamics ... "); tic
[~, XUNC] = ode15s(@(t, x) FofXU(x,   0    ), t, x0, opts); fprintf("completed in %2.2f seconds. \n", toc)
fprintf(" - Simulating LQR closed-loop dynamics ... "); tic
[~, XLQR] = ode15s(@(t, x) FofXU(x, uLQR(x)), t, x0, opts); fprintf("completed in %2.2f seconds. \n", toc)
fprintf(" - Simulating PPR closed-loop dynamics ... "); tic
[t, XPPR] = ode15s(@(t, x) FofXU(x, uPPR(x)), t, x0, opts); fprintf("completed in %2.2f seconds. \n", toc)
%% Animate solution
figure('Position', [311.6667 239.6667 1.0693e+03 573.3333]);
for i=1:length(t)
    
    Z = reshape(XUNC(i,:),nx,ny);
    
    subplot(2,2,1); grid on;
    [c,h] = contourf(X,Y,Z); clabel(c,h)
    xlabel('x, m'); ylabel('y, m'); axis equal
    
    subplot(2,2,2); surfc(X,Y,Z); zlim([-1.5 1.5])
    xlabel('x, m'); ylabel('y, m');
    
    Z = reshape(XPPR(i,:),nx,ny);
    
    subplot(2,2,3); grid on;
    [c,h] = contourf(X,Y,Z); clabel(c,h)
    xlabel('x, m'); ylabel('y, m'); axis equal
    
    subplot(2,2,4); surfc(X,Y,Z); zlim([-1.5 1.5])
    xlabel('x, m'); ylabel('y, m');
    
    drawnow
end

%% Plot solution
figure('Position', [311.6667 239.6667 1.0693e+03 573.3333]);
tlo = tiledlayout(3,4);ii=0;
for i=round(linspace(1,size(XUNC,1),4))
    ii=ii+1;h(ii) = nexttile;
    Z = reshape(XUNC(i,:),nx,ny);
    surfc(X,Y,Z); zlim([-1.5 1.5])
    xlabel('x, m'); ylabel('y, m');
end
for i=round(linspace(1,size(XLQR,1),4))
    ii=ii+1;h(ii) = nexttile;
    Z = reshape(XLQR(i,:),nx,ny);
    surfc(X,Y,Z); zlim([-1.5 1.5])
    xlabel('x, m'); ylabel('y, m');
end
for i=round(linspace(1,size(XPPR,1),4))
    ii=ii+1;h(ii) = nexttile;
    Z = reshape(XPPR(i,:),nx,ny);
    surfc(X,Y,Z); zlim([-1.5 1.5])
    xlabel('x, m'); ylabel('y, m');
end
set(h, 'Colormap', turbo, 'CLim', [-1 1])
cbh = colorbar(h(end));
cbh.Layout.Tile = 'east';
drawnow

end


function status = odeprog(t, y, flag)
    % ODEPROG Custom progress bar for ode solver
    % Use with odeset: opts = odeset('OutputFcn',@odeprog);
    persistent T1 nSteps lastPct lastUpdateTime
    global T0 

    status = false;

    switch flag
        case 'init'
            % Initialize progress bar
            elapsed = toc(T0);
            T1 = t(end);
            nSteps = 50; % Number of blocks in the progress bar
            lastPct = -1;
            lastUpdateTime = 0;
            fprintf(' |%s|  (elapsed: %5i s, remaining: ----- s)\n', repmat(' ',1,nSteps), round(elapsed));

        case ''
            % ODE solver step
            if isempty(t), return; end
            tNow = t(end);
            pct = min(100, max(0, 100 * tNow / T1));
            block = floor(pct / (100/nSteps));
            elapsed = toc(T0);
            eta = (elapsed / max(tNow,eps)) * (T1 - tNow); % avoid divide-by-zero
            needsUpdate = pct - lastPct >= 2 || block == nSteps;
            timeSinceLast = elapsed - lastUpdateTime;

            if needsUpdate || timeSinceLast >= 1
                bar = [repmat('-',1,block), repmat(' ',1,nSteps-block)];
                fprintf(repmat('\b',1,94));
                fprintf(' |%s|  (elapsed: %5i s, remaining: %5i s)\n', bar, round(elapsed), round(eta));
                if needsUpdate
                    lastPct = pct;
                end
                lastUpdateTime = elapsed;
            end

        case 'done'
            % Finalize
            % elapsed = toc(T0);
            % bar = repmat('-',1,nSteps);
            fprintf(repmat('\b',1,94));
            % fprintf(' |%s|  (elapsed: %5i s, remaining:     0 s)\n', bar, round(elapsed));
            clear T0 T1 nSteps lastPct lastUpdateTime
    end
end

