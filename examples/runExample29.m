function runExample29(numElements, r)
%runExample29 Runs the 2D Allen-Cahn FEM example with low-rank penalty
%
%   Usage:  runExample29(numElements)
%
%   Inputs:
%       numElements - number of finite elements in each direction
%       r           - reduced-order dimension
%
%   Description: Here we consider Allen-Cahn with Neumann BCs, insulated on 3
%   sides and subject to boundary control on the fourth. The PDE is
%
%     uₜ(x,y,t) = ε Δu(x,y,t) + u(x,y,t) - u(x,y,t)³
%     uᵧ(x,1,t) = u₃(t)  (Neumann control input on side CD)
%
%   where all sides are insulated except one subject to Neumann boundary
%   control. The PDE model can be discretized with the finite element method
%   [1,2], and the resulting finite-dimensional ODE model can be written as
%
%     M ẋ + K₁ x + K₃ (x⊗x⊗x) = B u,
%
%   for which we can compute a state feedback controller u(x) = K(x) using PPR.
%
%   This example demonstrates the benefits and importance of properly
%   exploiting sparsity for high-dimensional nonlinear control problems.
%   Sparsity is used in four key ways in this example:
%       1) Forming and storing the FEM model in generalized form
%       2) Computing the Riccati solutions using modern low-rank methods
%       3) Forming the projected ROM for the higher-order coefficients
%       4) Simulating the ODEs efficiently by using the sparse dynamics
%
%   By properly leveraging sparsity and applying dimensionality reduction,
%   the PPR method can be scaled to quite large systems, since it is
%   essentially a nonlinear update to existing powerful linear methods.
%   Since the first term is based on a Riccati equation, choosing the cost
%   function cleverly (to be low-rank) permits using the powerful modern
%   LR-ADI solvers from the M-M.E.S.S. package [3]. For the remaining
%   computations in terms of the Kronecker product, this example showcases
%   a new custom data structure that aids with for Kronecker
%   polynomials. The custom class, called sparseIJV, is essentially a sparse
%   array, but it stores only the indices, values, and size of the arrays.
%   For the level of sparsity (and array dimensions) that arise with
%   Kronecker polynomials, it allows major speedups. Lastly, during the ODE
%   simulation steps, evaluating the sparse dynamics, the Jacobian, and
%   leveraging the mass matrix all lead to major performance gains.
%
%   Combining all of these major performance considerations permits running
%   this model on a laptop in dimensions as high as n=103041 dimensions.
%   The PPR computation also runs in n=263169 dimension, but the 
%   closed-loop ode simulation runs out of memory. It can be run on the
%   workstation with 512GB RAM.
%
%                                           Total Script Time
%   +--------------+---------+----------------------+----------------------+---------------------+
%   | numElements  |    n    |   CPU Time Laptop    |   CPU Time Laptop    |   CPU Time Server   |
%   |              |         |     (16 GB RAM)      |     (32 GB RAM)      |     (512 GB RAM)    |
%   +--------------+---------+----------------------+----------------------+---------------------+
%   |      64      |   4225  |        36 sec        |         25 sec       |       27 sec        |
%   |     128      |  16641  |         4 min        |          2 min       |        4 min        |
%   |     256      |  66049  |        70 min        |         25 min       |       46 min        |
%   |     320      | 103041  |          OOM         |                      |                     |
%   |     512      | 263169  |          ---         |                      |                     |
%   +--------------+---------+----------------------+----------------------+---------------------+
%
%                                     PPR Control Computation Time
%   +--------------+---------+----------------------+----------------------+---------------------+
%   | numElements  |    n    |   CPU Time Laptop    |   CPU Time Laptop    |   CPU Time Server   |
%   |              |         |     (16 GB RAM)      |     (32 GB RAM)      |     (512 GB RAM)    |
%   +--------------+---------+----------------------+----------------------+---------------------+
%   |      64      |   4225  |        22 sec        |         14 sec       |       14 sec        |
%   |     128      |  16641  |        85 sec        |         48 sec       |       82 sec        |
%   |     256      |  66049  |         6 min        |          3 min       |        8 min        |
%   |     320      | 103041  |        11 min        |          5 min       |       11 min        |
%   |     512      | 263169  |         ---          |         11 min       |                     |
%   +--------------+---------+----------------------+----------------------+---------------------+
%
%   Reference: [1] D. M. Boskovic, M. Krstic, and W. Liu, "Boundary control
%              of an unstable heat equation via measurement of
%              domain-averaged temperature," IEEE Transactions on Automatic
%              Control, vol. 46, no. 12, pp. 2022–2028, 2001
%              [2] S. A. Ragab and H. E. Fayed, Introduction to finite
%              element analysis for engineers. Taylor & Francis Group, 2017
%              [3] J. Saak, M. Köhler, and P. Benner, "M-M.E.S.S. - the
%              matrix equation sparse solver library,"  v3.1, 2025.
%              doi: 10.5281/zenodo.632897. https://github.com/mpimd-csc/mmess
%
%   Part of the PPR repository.
%%
fprintf('Running Example 29\n')
if nargin < 2
    if nargin < 1
        numElements = 128;
    end
    r = 10;
end
animate = false;

%% Get dynamics
nx = numElements+1; ny = nx; n = nx*ny; m = 1;
fprintf(" Forming FEM model, n=%i ... ",n); tic
[E, f, g, ~, xyg] = getSystem29(numElements,.25,1,-1);
fprintf("completed in %2.2f seconds. \n", toc)
g{1} = g{1}(:,1);
boundaryLocs = xyg(g{1}>0, 1);
% g{1}(g{1}>0) = (-(2*boundaryLocs-1).^2+1) * max(g{1}); % Parabola
g{1}(g{1}>0) = sin(pi*boundaryLocs) * max(g{1}); % Parabola
G = @(x) g{1}; % insulate all sides, use control only on side CD

%% Compute controllers
% Setting the cost Q=C.'*C for LR-ADI
options.lrradi = true;
nc = 10; nds = round(linspace(1,n,nc));
C = sparse(1:nc,nds,1,nc,n); q = C.'*C;

% Get value function/controller
R = 0.1;
degree = 4;
options.C = C; options.E = E;
options.verbose = false; options.r = r;
fprintf(" Computing ppr() solution w/ lrradi, n=%i, r=%i, d=%i ... ",n,r,4); tic
[~, K] = ppr(f, g, q, R, degree, options);
fprintf(" completed in %2.2f seconds. \n", toc)

uPPR = @(x) kronPolyEval(K, x, degree-1);

%% Simulate closed-loop system
% Clear variables used in ode solve
clear sparseKronPolyEval T0; global T0;

% Set up ode function, mass matrix, and Jacobian
% (to maximize efficient sparsity usage)
% FofXU = @(x,u) kronPolyEval(f,x) + g{1} * u;      % normal kronPolyEval (about 10% slower than custom one)
FofXU = @(x,u) sparseKronPolyEval(f,x) + g{1} * u;  % custom kronPolyEval with additional slight improvement
opts_openloop = odeset(Mass=E, Jacobian=f{1},           OutputFcn=@odeprog);

opts_closloop = odeset(Mass=E, Jacobian=f{1}+g{1}*sparse(K{1}), OutputFcn=@odeprog);

% Set up grid and annulus initial condition
X = reshape(xyg(:,1),nx,ny); Y = reshape(xyg(:,2),nx,ny);

% Initial Condition 1: Sine wave
% x0 = .25*(sin(4*pi*X) + cos(3*pi*Y)) + .1; x0 = x0(:);

% Initial Condition 2: A thin annulus (circle)
% Rad = sqrt((X - 0.5).^2 + (Y - 0.5).^2); % get radius values for grid points
% D = 0.75; thickness = 0.03; % annulus with outer diameter 0.75 and thickness 0.015
% x0 = (Rad <= D/2 & Rad >= D/2-thickness );
% x0 = 0.5*x0(:)+0.5; % x0 is made as a logical by the last line

% Initial Condition 3: Stanford Bunny
x0 = getStanfordBunnyIC(X, Y);
x0 = x0(:)*0.5+0.5;

tmax = 5; t = (0:0.005:1).^3 * tmax; % specify for consistent plotting

t = [0, 0.005, 0.25, tmax];

% Run and time simulations
fprintf(" - Simulating open-loop dynamics ... ");
T0 = tic;
[~, XUNC] = ode15s(@(t, x) FofXU(x,   0    ), t, x0, opts_openloop);
fprintf("completed in %2.2f seconds. \n", toc(T0))
fprintf(" - Simulating PPR closed-loop dynamics ... ");
T0 = tic;
[~, XPPR] = ode15s(@(t, x) FofXU(x, uPPR(x)), t, x0, opts_closloop);
fprintf("completed in %2.2f seconds. \n", toc(T0))

%% Compute integrated costs
% costUNC = trapz(t, sum((XUNC.^2).*diag(C.'*C).', 2));
% for i=1:length(t); UxPPR(i,1) = uPPR(XPPR(i,:).'); end
% costPPR = trapz(t, sum((XPPR.^2).*diag(C.'*C).', 2) + R*UxPPR.^2);

%% Plot solution
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
[~, idx1] = min(abs(t - .005));
[~, idx2] = min(abs(t - .25));
figure('Position', [312 600 864*1.25 195*1.25]);
tiledlayout(1,4);ii=0;
for i=[1, idx1, idx2, length(t)]
    ii=ii+1;h(ii) = nexttile;
    Z = reshape(XUNC(i,:),nx,ny);
    surfc(X,Y,Z,'EdgeAlpha',0); zlim([-1.5 1.5])
    xlabel('$z_1$','Interpreter','latex'); ylabel('$z_2$','Interpreter','latex');
    title(sprintf('$t=%4.3f$',t(i)),'Interpreter','latex')
end
drawnow
exportgraphics(gcf,'plots/example29_UNC.png', 'ContentType', 'image')

figure('Position', [312 300 864*1.25 195*1.25]);
tiledlayout(1,4);ii=0;
for i=[1, idx1, idx2, length(t)]
    ii=ii+1;h(ii) = nexttile;
    Z = reshape(XPPR(i,:),nx,ny);
    surfc(X,Y,Z,'EdgeAlpha',0); zlim([-1.5 1.5])
    xlabel('$z_1$','Interpreter','latex'); ylabel('$z_2$','Interpreter','latex');
    title(sprintf('$t=%4.3f$',t(i)),'Interpreter','latex')
end
drawnow
exportgraphics(gcf,'plots/example29_PPR.png', 'ContentType', 'image')


%% Animate solution
if animate
    figure('Position', [311.6667 239.6667 1.0693e+03 573.3333]);
    for i=1:5:length(t)
        
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
end

clear sparseKronPolyEval T0
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
            fprintf(' |%s|  (elapsed: %5i s, remaining: %5i s)', bar, round(elapsed), min(round(eta),99999));
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
if isempty(n) || n ~= length(z)
    n = length(z);
    [F3i, F3j, F3v] = find(f{3});
    [I1, I2, I3] = ind2sub([n n n], F3j);
    I1 = uint32(I1); I2 = uint32(I2); I3 = uint32(I3);
end

% Evaluate linear and quadratic terms normally
x = f{1}*z;
% x = x + f{2}*kron(z,z); % commented out because f{2}=0

% Efficient sparse evaluation of f{3}*(z⊗z⊗z)
zprod = z(I1) .* z(I2) .* z(I3); % uint necessary for large models, double has rounding error
x = x + accumarray(F3i, F3v .* zprod, size(x));


end


function Zq = getStanfordBunnyIC(Xq, Yq)
%getStanfordBunnyIC Creates a 2D initial condition resembling the Stanford
%bunny for a square domain

x = [-77; -77; -77; -76; -75; -74; -71; -71; -72; -72; -73; -73; -71; -71; -69; -68; -67; -65; -62; -60; -58; -56; -51; -50; -45; -40; -44; -50; -54; -56; -53; -48; -42; -41; -38; -37; -35; -31; -29; -25; -24; -23; -22; -21; -20; -19; -17; -16; -10; -8; -4; 0; 7; 13; 18; 21; 22; 23; 25; 27; 28; 32; 37; 40; 42; 45; 52; 55; 61; 62; 69; 71; 73; 75; 77; 77; 78; 77; 75; 72; 70; 67; 64; 62; 63; 62; 62; 61; 60; 58; 58; 57; 54; 52; 50; 46; 41; 39; 35; 32; 28; 26; 20; 19; 12; 11; 7; 4; -1; -7; -14; -18; -20; -24; -27; -31; -32; -31; -30; -28; -23; -17; -14; -8; -2; 3; 4; 6; 7; 8; 8; 6; 6; 3; 2; -1; -3; -5; -11; -18; -20; -26; -29; -34; -44; -45; -46; -49; -50; -54; -57; -60; -60; -62; -62; -61; -66; -68; -71; -72; -74; -74; -74; -74; -74; -76; -77];
y = [16; 11; 9; 7; 5; 4; -1; -2; -7; -12; -16; -20; -25; -26; -30; -33; -34; -37; -40; -42; -44; -46; -48; -49; -51; -58; -63; -66; -70; -74; -76; -77; -77; -77; -77; -77; -77; -77; -77; -76; -76; -76; -76; -76; -76; -76; -76; -76; -76; -75; -75; -77; -77; -76; -75; -76; -76; -76; -76; -76; -76; -76; -76; -75; -75; -75; -73; -71; -70; -69; -65; -63; -62; -57; -51; -50; -44; -43; -39; -38; -37; -36; -34; -29; -26; -19; -17; -12; -11; -8; -4; -3; 2; 4; 8; 12; 15; 16; 18; 20; 21; 21; 21; 21; 21; 21; 19; 19; 18; 17; 16; 17; 17; 18; 20; 24; 27; 34; 36; 38; 39; 43; 43; 46; 49; 52; 52; 56; 58; 62; 64; 69; 70; 73; 76; 76; 76; 76; 72; 67; 66; 61; 57; 55; 60; 64; 66; 69; 70; 70; 70; 67; 66; 60; 55; 52; 44; 44; 42; 41; 36; 33; 32; 29; 25; 20; 16];

lim = 100;
[X, Y] = meshgrid(-lim:lim, -lim:lim);

Z = X*0;

for i=1:numel(Z)
    if inpolygon(X(i),Y(i),x, y)
        Z(i) = 1000;
        for j=1:length(x)
            if norm([X(i)-x(j), Y(i)-y(j)]) < Z(i)
                Z(i) = norm([X(i)-x(j), Y(i)-y(j)]);
            end
        end
        if Z(i) > 10
            Z(i) = 0;
        else
            Z(i) = 1;
        end
    else
        Z(i) = 0;
    end
end

% surf(X, Y, Z);
Zq = interp2(X./lim./2+.5,Y./lim./2+.5,Z,Xq,Yq);
end

