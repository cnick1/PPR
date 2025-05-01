function runExample29_LR(numElements, r)
%runExample29_LR Runs the 2D unsteady nonlinear heat equation FEM control example.
%
%   Usage:  runExample29_LR(numElements)
%
%   Inputs:
%       numElements - number of finite elements in each direction
%       r           - reduced-order dimension
%
%   Description: The model, inspired by [1], describes an unstable heat
%   equation problem modeling heat generated in a resistive electrical
%   material. The resulting model takes the form
%
%     uₜ(x,y,t) = uₓₓ(x,y,t) + uᵧᵧ(x,y,t) + λ u(x,y,t)
%     uᵧ(x,1,t) = u₃(t)  (Neumann control input on side CD)
%
%   where all sides are insulated except one subject to Neumann boundary
%   control. The FEM model can be written as
%
%     M ẋ + K₁ x + K₃ (x⊗x⊗x) = B u,
%     y = C x.
%   
%   for which we can compute a controller u(x) = K(x) using PPR.
%
%   This example demonstrates the benefits and importance of properly 
%   exploiting sparsity. Sparsity is used in four key ways in this example:
%       1) Forming and storing the FEM model in generalized form
%       2) Computing the Riccati solutions using low-rank methods
%       3) Forming the projected ROM for the higher-order coefficients
%       4) Simulating the odes by using the sparse dynamics
%
%   This example illustrates how scalable the PPR method is, since it is
%   essentially a nonlinear update to existing powerful linear methods.
%   Since the first term is based on a Riccati equation, choosing the cost
%   function cleverly (to be low-rank) permits using the powerful modern
%   LR-RADI solvers from the M-M.E.S.S. package. For the remaining
%   computations in terms of the Kronecker product, this example showcases
%   a new custom data structure that is necessary for Kronecker
%   polynomials. The custom class, called sparseIJV, is essentially a sparse 
%   array, but it stores only the indices, values, and size of the arrays.
%   For the level of sparsity (and array dimensions) that arise with
%   Kronecker polynomials, it allows major speedups. Lastly, during the ode 
%   simulation steps, evaluating the sparse dynamics, the Jacobian, and 
%   leveraging the mass matrix all lead to major performance gains.
%
%   Combining all of these major performance considerations permits running
%   this model on a laptop in dimensions as high as n=16641 dimensions, and
%   on a workstation with 512 GB RAM up to n=66049 dimensions. Here is a
%   summary of the performance on a laptop vs server for different dimensions:

%                            Total Script Time
%   +--------------+---------+----------------------+---------------------+
%   | numElements  |    n    | CPU Time Laptop      | CPU Time Server     |
%   |              |         | (16 GB RAM)          | (512 GB RAM)        |
%   +--------------+---------+----------------------+---------------------+
%   |      64      |  4225   |        60 sec        |       50 sec        |
%   |     128      | 16641   |        30 min        |       10 min        |
%   |     200      | 40401   |          --          |        2  hr        |
%   |     256      | 66049   |          --          |        6  hr        |
%   +--------------+---------+----------------------+---------------------+
%
%                        PPR Control Computation Time
%   +--------------+---------+----------------------+---------------------+
%   | numElements  |    n    | CPU Time Laptop      | CPU Time Server     |
%   |              |         | (16 GB RAM)          | (512 GB RAM)        |
%   +--------------+---------+----------------------+---------------------+
%   |      64      |  4225   |        20 sec        |       15 sec        |
%   |     128      | 16641   |        90 sec        |       60 sec        |
%   |     200      | 40401   |        40 min        |        5 min        |
%   |     256      | 66049   |          --          |       10 min        |
%   +--------------+---------+----------------------+---------------------+
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
if nargin < 2
    if nargin < 1
        numElements = 50;          
    end
    r = 10;
end

%% Get dynamics
nx = numElements+1; ny = nx; n = nx*ny; m = 1;
fprintf(" Forming FEM model, n=%i ... ",n); tic
[E, f, g, ~, xyg] = getSystem29(numElements,.75,1,-1);
fprintf("completed in %2.2f seconds. \n", toc)
g{1} = g{1}(:,3); G = @(x) g{1}; % insulate all sides, use control only on side CD

%% Compute controllers
% Setting the cost Q=C.'*C for LR-ADI
options.lrradi = 1;
nc = 10; nds = round(linspace(1,n,nc));
C = sparse(1:nc,nds,sqrt(0.1),nc,n); q = C.'*C;

% Get value function/controller
R = 1; degree = 4;
options.C = C; options.E = E;
options.verbose = false; options.r = r;
fprintf(" Computing ppr() solution w/ lrradi, n=%i, r=%i, d=%i ... ",n,r,4); tic
[~, K] = ppr(f, g, q, R, degree, options);
fprintf(" completed in %2.2f seconds. \n", toc)

uLQR = @(x) kronPolyEval(K, x, 1);
uPPR = @(x) kronPolyEval(K, x, degree-1);

%% Simulate closed-loop system
% Clear variables used in ode solve
clear F3i F3j F3v I1 I2 I3 T0; global T0;

% Set up ode function, mass matrix, and Jacobian
% (to maximize efficient sparsity usage)
% FofXU = @(x,u) kronPolyEval(f,x) + g{1} * u;      % normal kronPolyEval (about 10% slower than custom one)
FofXU = @(x,u) sparseKronPolyEval(f,x) + g{1} * u;  % custom kronPolyEval with additional slight improvement
opts_openloop = odeset(Mass=E, Jacobian=f{1},           OutputFcn=@odeprog); % option 3
opts_closloop = odeset(Mass=E, Jacobian=f{1}+g{1}*K{1}, OutputFcn=@odeprog); % option 3

% Set up grid and annulus initial condition
X = reshape(xyg(:,1),nx,ny); Y = reshape(xyg(:,2),nx,ny);
% x0 = .25*(sin(4*pi*X) + cos(3*pi*Y)) + .1; x0 = x0(:);
R = sqrt((X - 0.5).^2 + (Y - 0.5).^2); % get radius values for grid points
x0 = (R >= 0.360 & R <= 0.375); % annulus with diameter 0.75 and thickness 0.015
x0 = 0.5*x0(:); % x0 is made as a logical by the last line

tmax = 5; t = (0:0.005:1).^3 * tmax; % specify for consistent plotting

% Run and time simulations
fprintf(" - Simulating open-loop dynamics ... ");       T0 = tic;
[~, XUNC] = ode15s(@(t, x) FofXU(x,   0    ), t, x0, opts_openloop); fprintf("completed in %2.2f seconds. \n", toc(T0))
fprintf(" - Simulating LQR closed-loop dynamics ... "); T0 = tic;
[~, XLQR] = ode15s(@(t, x) FofXU(x, uLQR(x)), t, x0, opts_closloop); fprintf("completed in %2.2f seconds. \n", toc(T0))
fprintf(" - Simulating PPR closed-loop dynamics ... "); T0 = tic;
[~, XPPR] = ode15s(@(t, x) FofXU(x, uPPR(x)), t, x0, opts_closloop); fprintf("completed in %2.2f seconds. \n", toc(T0))

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
    title(sprintf('t=%7.6f',t(i)))
    drawnow
end

clear n F3i F3j F3v I1 I2 I3 T0
% whos f
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
        fprintf(' |%s|  (elapsed: %5i s, remaining: ----- s)', repmat(' ',1,nSteps), round(elapsed));
        
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
            fprintf(repmat('\b',1,93));
            fprintf(' |%s|  (elapsed: %5i s, remaining: %5i s)', bar, round(elapsed), round(eta));
            if needsUpdate
                lastPct = pct;
            end
            lastUpdateTime = elapsed;
        end
        
    case 'done'
        % Finalize
        % elapsed = toc(T0);
        % bar = repmat('-',1,nSteps);
        fprintf(repmat('\b',1,93));
        % fprintf(' |%s|  (elapsed: %5i s, remaining:     0 s)\n', bar, round(elapsed));
        clear T0 T1 nSteps lastPct lastUpdateTime
end
end

function [x] = sparseKronPolyEval(f,z)
%sparseKronPolyEval Evaluate a Kronecker polynomial with sparse optimization
% Note: This is ONLY for runExample29:
%   - Assumes f{2} is zero, f{1} and f{3} are only other coefficients
%   - Assumes f{3} is sparse and avoids forming kron(z,z,z)
%
% Inputs:   f - coefficients cell array, f{1}, f{2}, f{3}
%           z - value to calculate the polynomial at (vector)
%
% Output:   x = f{1}*z + f{3}*(z⊗z⊗z)
%%

% Use persistent variables to only compute sparsity pattern once
% (major speed-up)
persistent n F3i F3j F3v I1 I2 I3
if isempty(n)
    n = length(z);
    [F3i, F3j, F3v] = find(f{3});
    [I1, I2, I3] = ind2sub([n n n], F3j);
end

% Evaluate linear and quadratic terms normally
x = f{1}*z;
% x = x + f{2}*kron(z,z); % commented out because f{2}=0

% Efficient sparse evaluation of f{3}*(z⊗z⊗z)
zprod = z(I1) .* z(I2) .* z(I3);
x = x + accumarray(F3i, F3v .* zprod, size(x));

end

