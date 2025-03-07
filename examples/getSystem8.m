function [f, g, h] = getSystem8(numElements, elementOrder)
%getSystem8  Generates a cubic finite element heat equation model system 
%   The system is a finite element model for a nonlinear heat equation, i.e.
%   a reaction-diffusion equation. The function returns a finite element
%   model with either linear or quadratic elements (elementOrder 1 or 2). 
%   The state space dimension is n=numElements-1 for linear elements or
%   n=2*numElements-1 for quadratic elements, due to the fixed end boundary
%   conditions.
%
%   Usage:   [f,g,h] = getSystem8()
%         or [f,g,h] = getSystem8(numElements,elementOrder)
%
%   Inputs:
%       numElements    - number of elements to discretize the domain with
%                        (default = 4)
%       elementOrder   - select either:
%                          • 1 = linear elements (default)
%                          • 2 = quadratic elements (not fully implemented)
%
%   Outputs:    f,g,h  - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%
%   Description: The eaction-diffusion model used in [1] from [2-4] is a
%   represented by the nonlinear heat equation PDE 
%
%     uₜ(x,t) = uₓₓ(x,t) + uₓ(x,t) + 1/8 u(x,t) + u(x,t)³
%     u(0,t) = 0                     (Dirichlet BC)
%     u(1,t) = 0                     (Dirichlet BC)
% 
%   which we augment with a control input. After FEM discretization, the 
%   finite element equations can be written as
%
%     M ẋ + K(x) x = B(x) u,
%     y = C x.
%
%   The terms K(x) and B(x) can be approximated with Taylor series
%   expansions, which leads to the Kronecker product representation
%
%     M ẋ + K₁ x + K₂ (x ⊗ x) + K₃ (x ⊗ x ⊗ x) + ...
%       = B₀ u + B₁ (x ⊗ u) + B₂ (x ⊗ x ⊗ u) + ...,
%     y = C x.
%
%   Reference: [1] N. A. Corbin and B. Kramer, “Scalable computation of
%               𝓗∞ energy functions for polynomial drift nonlinear
%               systems,” in 2024 American Control Conference (ACC), Jul.
%               2024, pp. 2506–2511. doi: 10.23919/acc60939.2024.10644363
%              [2] M. Embree, “Unstable modes in projection-based
%               reduced-order models: how many can there be, and what do they
%               tell you?,” Systems & Control Letters, vol. 124, pp. 49–59,
%               Feb. 2019, doi: 10.1016/j.sysconle.2018.11.010
%              [3] J. Galkowski, “Nonlinear instability in a semiclassical
%               problem,” Communications in Mathematical Physics, vol. 316,
%               no. 3, pp. 705–722, Oct. 2012, doi: 10.1007/s00220-012-1598-5
%              [4] B. Sandstede and A. Scheel, “Basin boundaries and
%               bifurcations near convective instabilities: a case study,”
%               Journal of Differential Equations, vol. 208, no. 1, pp.
%               176–193, Jan. 2005, doi: 10.1016/j.jde.2004.02.016
%
%   Part of the PPR repository.
%%
vec = @(X) X(:);
if nargin < 2
    elementOrder = 1;
    if nargin < 1
        numElements = 4;
    end
end

if elementOrder == 2; error("Not fully implemented!"); end

%% Prepare linear domain mesh parameters
nel = numElements;            % number of elements
nng = elementOrder * nel + 1; % number of nodes
nvpn = 1;                     % number of variables per node
nnpe = elementOrder + 1;      % number of nodes per element
nvpe = nnpe * nvpn;           % number of variables per element
nvg = nng * nvpn;             % number of variables in global mesh

%% Specify material and geometry properties
% Define geometry and material properties
% I am using notation for the model problem (5.1) in Ragab [2]
L = 30;                       % wire length
gamma = 1;                    % mass property
p = 1;                        % nondimensionalized diffusion coefficient
q = -1/8;                     % -λ, to do with heat generation due to electrical resistance

%% Generate linear "mesh"
xg = linspace(0, L, nng);     % node locations
he = xg(nnpe) - xg(1);        % element length

%% Assemble linear global matrices (mass and stiffness)
% Define mass matrix for one element
if elementOrder == 1
    Me = he / 6 * ...
        [2, 1;
        1, 2];
    
    % Define stiffness matrix for one element
    K1e = 1 / he * [1, -1;
        -1, 1] ...
        + q * he / 6 * [2, 1;
        1, 2] ...
        +1/2 * [-1 1;
        -1 1];
elseif elementOrder == 2
    Me = he / 30 * ...
        [4, 2, -1;
        2, 16, 2;
        -1, 2, 4];
    
    % Define stiffness matrix for one element
    K1e = 1 / (3 * he) * [7, -8, 1;
        -8, 16, -8;
        1, -8, 7] ...
        - he / 240 * [4, 2, -1;
        2, 16, 2;
        -1, 2, 4];
end
% Initialize and stack/assemble global matrix
Mg = sparse(nvg, nvg);
K1g = sparse(nvg, nvg);

if elementOrder == 1
    for ii = 1:nel
        Mg(ii:(ii + nvpe - 1), ii:(ii + nvpe - 1)) = Mg(ii:(ii + nvpe - 1), ii:(ii + nvpe - 1)) + Me;
        K1g(ii:(ii + nvpe - 1), ii:(ii + nvpe - 1)) = K1g(ii:(ii + nvpe - 1), ii:(ii + nvpe - 1)) + K1e;
    end
elseif elementOrder == 2
    for i = 1:nel
        ii = 2 * i - 1;
        Mg(ii:(ii + nvpe - 1), ii:(ii + nvpe - 1)) = Mg(ii:(ii + nvpe - 1), ii:(ii + nvpe - 1)) + Me;
        K1g(ii:(ii + nvpe - 1), ii:(ii + nvpe - 1)) = K1g(ii:(ii + nvpe - 1), ii:(ii + nvpe - 1)) + K1e;
    end
end

%% Assemble quadratic global matrix
% Form all zero quadratic global matrix directly
K2g = sparse(nvg, nvg ^ 2);

%% Assemble cubic global matrix
% Define stiffness matrix for one element
if elementOrder == 1
    K3e = -he / 20 * [4, 3, 0, 2, 0, 0, 0, 1;
        1, 2, 0, 3, 0, 0, 0, 4];
elseif elementOrder == 2
    K3e = -he / 1260 * ...
        [92, 96, -21, 0, 96, -24, 0, 0, 6, 0, 0, 0, 0, 32, -48, 0, 0, -12, 0, 0, 0, 0, 0, 0, 0, 0, -7;
        32, 96, -12, 0, 96, -96, 0, 0, -12, 0, 0, 0, 0, 512, 96, 0, 0, 96, 0, 0, 0, 0, 0, 0, 0, 0, 32;
        -7, -12, 6, 0, -48, -24, 0, 0, -21, 0, 0, 0, 0, 32, 96, 0, 0, 96, 0, 0, 0, 0, 0, 0, 0, 0, 92];
end

% Initialize and stack/assemble global matrix
K3g = sparse(nvg, nvg ^ 3);

if elementOrder == 1
    for i = 0:nel - 1 % start from zero so you don't have to subtract 1 every time
        ii = i * nvpn + 1;
        idxs = (nvg ^ 2 + nvg + 1) * i * nvpn ... % Starting index shift depending on element iteration
            + vec(( ...
            vec(([1:nvpe] + [0:nvg:nvg * (nvpe - 1)]')')' ... % (basically the quadratic indices)
            + [0:nvg ^ 2:nvg ^ 2 * (nvpe - 1)]' ... % Add secondary skips into sequence (add row to column and then vec)
            )')';
        
        % "stack" element matrices into global matrix
        K3g(ii:(ii + nvpe - 1), idxs) = K3g(ii:(ii + nvpe - 1), idxs) + K3e;
    end
elseif elementOrder == 2
    for i = 1:nel
        ii = 2 * i - 1;
        idxs = (nvg ^ 2 + nvg + 1) * (i - 1) * nvpn * 2 ... % Starting index shift depending on element iteration
            + vec(( ...
            vec(([1:nvpe] + [0:nvg:nvg * (nvpe - 1)]')')' ... % (basically the quadratic indices)
            + [0:nvg ^ 2:nvg ^ 2 * (nvpe - 1)]' ... % Add secondary skips into sequence (add row to column and then vec)
            )')';
        
        % "stack" element matrices into global matrix
        K3g(ii:(ii + nvpe - 1), idxs) = K3g(ii:(ii + nvpe - 1), idxs) + K3e;
    end
end

%% RHS
% TODO
% RB0 = sparse(TotalDOFs, 2);
%         RB1 = sparse(TotalDOFs, 2 * TotalDOFs);
%         RB2 = sparse(TotalDOFs, 2 * TotalDOFs ^ 2);
%         RB3 = sparse(TotalDOFs, 2 * TotalDOFs ^ 3);

%         RB0(TotalDOFs - 2, 2) = 1; % Force in x direction
%         RB0(TotalDOFs - 1, 1) = 1; % Force in y direction
%         RB0(TotalDOFs, :) = 0; % Moment in z direction

RB0 = generate_B_matrix(4, (nvg - 1) / 4);

%% Impose boundary conditions
fixedDOFs = [1, nvg]; % First and last nodes fixed
freeDOFs = setdiff(1:nvg, fixedDOFs);

% Reduced system method
K1g = K1g(freeDOFs, freeDOFs);
Mg = Mg(freeDOFs, freeDOFs);
RB0 = RB0(freeDOFs, :);

% K2G could clean up
fixedDOFsSquared = vec(fixedDOFs.' + (0:nvg:nvg ^ 2 - 1));
fixedDOFsSquared = unique([fixedDOFsSquared; vec((1:nvg).' + (fixedDOFs - 1) * nvg)]);

freeDOFsSquared = setdiff(1:nvg ^ 2, fixedDOFsSquared);
K2g = K2g(freeDOFs, freeDOFsSquared);

% K3G could clean up
fixedDOFsCubed = vec(vec((fixedDOFs - 1) * nvg + [1:nvg].') + (0:nvg ^ 2:nvg ^ 3 - 1));
fixedDOFsCubed = [fixedDOFsCubed; vec((1:nvg ^ 2).' + (fixedDOFs - 1) * nvg ^ 2)]; % Top rows zero
fixedDOFsCubed = [fixedDOFsCubed; vec(vec(fixedDOFs.' + (0:nvg:nvg ^ 2 - 1)) + (0:nvg ^ 2:nvg ^ 3 - 1))];
fixedDOFsCubed = sort(unique(fixedDOFsCubed));

freeDOFsCubed = setdiff(1:nvg ^ 3, fixedDOFsCubed);
K3g = K3g(freeDOFs, freeDOFsCubed);

%% Convert to state-space representation
n = length(Mg);

McholL = chol(Mg).'; % Use Cholesky factor for inverting rather than inv()

F1 = -McholL.' \ (McholL \ K1g);

G0 = McholL.' \ (McholL \ RB0);

C = RB0.';

% Construct N₂
F2 = -McholL.' \ (McholL \ K2g);

% Construct N₃
F3 = -McholL.' \ (McholL \ K3g);

%% Format outputs
f = {full(F1), F2, F3};
g = {full(G0)};
h = {full(C)};

end

function B = generate_B_matrix(m, i)
if rem(i, 1) ~= 0
    error("Incorrect number of elements/inputs")
end
n = i * m + 1;
B = sparse(n, m);
for j = 0:m - 1
    B(i * j + 1:i * j + 1 + i, j + 1) = 1;
end

B = B ./ (i + 1);
end
