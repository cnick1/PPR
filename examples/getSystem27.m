function [f, g, h, xg] = getSystem27(numElements, lambda)
%getSystem27 Generates a finite element model for a 1D unstable heat equation for testing controllers.
%   The PDE described in [1] models a thin wire with electricity running
%   through it. The resistance of the wire causes internal heat generation,
%   which is counteracted by the convection/radiation of heat from the wire.
%
%   Usage:   [f,g,h] = getSystem27()
%
%   Inputs:
%       numElements    - number of elements to discretize the domain with
%                        (default = 19)
%       lambda         - coefficient on the reaction term; drives
%                        instability for lambda > 0
%
%   Outputs:    f,g,h  - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%                  xg  - global node locations
%
%   Description: The wire is modeled by heat conduction in a thin rod of
%   small constant cross section. All of the material properties are taken
%   as constant EXCEPT the electrical resistivity, which varies linearly
%   with respect to the temperature of the wire. The resulting LINEAR
%   nondimensionalized model in [1] is of the form
%
%     uₜ(x,t) = uₓₓ(x,t) + λ u(x,t)
%     uₓ(0,t) = 0                     (Neumann BC)
%     uₓ(1,t) = u(t)                  (Neumann boundary control input)
%
%   hence the model is a simple reaction-diffusion heat equation model.
%   One boundary features a Nuemann BC (insulated end), whereas the other
%   boundary features a Neumann control input. (Nuemann control is easy to
%   implement since it is just in the secondary variable). After finite
%   element discretization, the finite element equations for the reaction-
%   diffusion problem can be written as
%
%     M ẋ + K x = B u,
%     y = C x.
%
%   which can of course be put in the standard state-space form by
%   multiplying by the inverse of the mass matrix. In [1], the main results
%   correspond to Dirichlet boundary control, so here we made some changes
%   to more easily implement it with Neumann control (mentioned in [1] but
%   not demonstrated).
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
if nargin < 2
    lambda = 0;
    if nargin < 1
        numElements = 20;
    end
end

%% Prepare linear domain mesh parameters
nel = numElements;         % number of elements
nng = nel + 1;             % number of nodes global
nvpn = 1;                  % number of variables per node
nnpe = 2;                  % number of nodes per element
nvpe = nnpe * nvpn;        % number of variables per element
nvg = nng * nvpn;          % number of variables in global mesh

%% Specify material and geometry properties
% Define geometry and material properties
% I am using notation for the model problem (5.1) in Ragab [2]
L = 1;                     % wire length
gamma = 1;                 % mass property
p = 1;                     % nondimensionalized diffusion coefficient
q = -lambda;               % -λ, to do with heat generation due to electrical resistance

%% Generate linear "mesh"
xg = linspace(0, L, nng);  % node locations
he = xg(2) - xg(1);        % element length

%% Generate element matrices Me, Ke, and Re and Assemble global system Mg, Kg, Rg
% Define element matrix components
L0 = he/6 * [2, 1;
    1, 2];

L1 = 1/he * [1, -1;
    -1,  1];

% Define element mass and stiffness matrices
Me = gamma*L0;
Ke = p*L1 + q*L0;

% Initialize and stack/assemble global matrix
Mg = zeros(nvg, nvg);
Kg = zeros(nvg, nvg);
Fg = zeros(nvg, 1);
Rg = zeros(nvg, 1);
for ie = 1:nel
    nodes = ie:(ie + nvpe - 1);
    
    % Assemble element matrices Me, Ke into global matrix Mg, Kg
    Mg(nodes,nodes) = Mg(nodes, nodes) + Me;
    Kg(nodes,nodes) = Kg(nodes, nodes) + Ke;
end

%% Prescribe and apply boundary conditions
% x=0 end is insultated; corresponding node is #1
% Rg(1) = Rg(1) + 0; % Rg is already zero

% x=L end has Neumann boundary control
Rg(end) = Rg(end) + 1; % 1 is just a placeholder for B; signal will be u(t)

%% Put into standard control system form
% The current ODE we have for the global system is Mg u̇ + Kg u = Fg + Rg.
% We wish to put this in the standard state-space form for control systems:
%       ẋ = A x + B u
%       y = C x
% Multiplying my M⁻¹ achieves this. The new system will be
%       ̇u = -Mg⁻¹Kg u + Mg⁻¹Fg + Mg⁻¹Rg
%        := A u  + F
% For computational purposes, we will do this inversion by computing the
% Cholesky factor of M; if we were to do M\K, M\F, and M\R these would all
% solve using the Cholesky factor, so we are recycling it. (We could
% simulate faster using sparse matrices and specifying mass matrix as an
% ode option; the inversion destroys the sparsity.)

Mchol = chol(Mg).'; % Use Cholesky factor for inverting
A = -Mchol.' \ (Mchol \ Kg);
% F =  Mchol.' \ (Mchol \ (Fg + Rg));

% Now, changing notation a bit, we want to put the system in the form
%       ̇x = f(x) + g(x) u
% so the old u becomes x and we introduce this new u for controls. f(x) is
% easy, it is just f(x) = A x. g(x) is going to be a B matrix; it is going
% to be computed by taking the F matrix; the control input u is really
% going through the qs variables in Rg, so if you add "forcing" fgen later
% this would have to change.
B = Rg;
B =  Mchol.' \ (Mchol \ B);

%% Placeholder for nonlinear terms
F2 = sparse(nvg,nvg^2);
F3 = sparse(nvg,nvg^3);

%% Format outputs
f = {A, F2, F3};
g = {B};
h = {1/nvg * ones(1,nvg)}; % average over domain


end

