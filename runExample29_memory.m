function runExample29_memory(numElements)
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
if nargin < 1
    numElements = 32;
end
% numElements = 512;
r = 10;

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
options.verbose = false; options.reducedDimension = r;
fprintf(" Computing ppr() solution w/ lrradi, n=%i, r=%i, d=%i ... ",n,r,4); tic
[v, K] = ppr(f, g, q, R, degree, options);
fprintf(" completed in %2.2f seconds. \n", toc)

V2 = v{2};
% v
% whos V2 v

fprintf('The full-order Riccati solution requires %f GB of RAM;\n the low-rank solution requires %f GB of RAM.\n',n^2 * 8/1e9,nnz(V2.Z)*8/1e9)
fprintf('The full-order v4 solution requires %f GB of RAM;\n the reduced v4 solution requires %f GB of RAM.\n',n^4 * 8/1e9,(nnz(v.ReducedValueCoefficients{4}) + nnz(v.Tinv))*8/1e9)

end
