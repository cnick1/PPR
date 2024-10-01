function [v,K,options] = ppr(f, g, q, r, degree, options)
%ppr  Compute a polynomial approximation to the value function for a polynomial
% control-affine dynamical system.
%
%   Usage: v = ppr(f, g, q, R, degree)
%
%       H∞ balancing energy functions can be computed as
%           [v] = ppr(f, g, cellfun(@(x) x * (-eta), h2q(h), 'un', 0), -1, degree);
%           [w] = ppr(f, g, h2q(h), 1/eta, degree);
%
%   Inputs:
%       f,g     - cell arrays containing the polynomial coefficients
%                 for the drift and input.
%       q       - cell arrays containing the polynomial coefficients for the
%                 state penalty in the cost. At minimum q{2} = Q(:) where Q is
%                 the quadratic weighting, i.e. x.'*Q*x (can optionally pass Q
%                 matrix directly in this case)
%       r       - cell arrays containing the polynomial coefficients for the
%                 control penalty in the cost. At minimum r{1} = R where R is
%                 the quadratic weighting, i.e. u.'*R*u (can optionally pass R
%                 matrix directly in this case)
%       degree  - desired degree of the computed value function. A degree d
%                 energy function uses information from f,g,q up-to degree d-1.
%                 The default choice of d is lf+1, where lf is the degree of
%                 the drift.
%       options - struct containing additional options:
%           - verbose: print runtime information
%         * The next set of options are used for accelerating the computations.
%           We can use linear balanced truncation to compute a reduced-order model
%           of dimension r, and then use the reduced dynamics to compute the energy
%           functions in r dimensions rather than n; this can present a huge
%           computational speedup.
%           - r: dimension "r" of reduced order model. (Defaults to n)
%           - fr,gr,qr: reduced dynamics; if for example the past energy
%             function has already been computed, instead of recomputing the
%             reduced dynamics, they can be passed in directly.
%           - eta: H-infinity balancing parameter. Defaults to open-loop balacing
%                 (eta = 0).
%
%
%   Output:
%       v       - cell array containing the polynomial value function coefficients
%       K       - cell array containing the polynomial (sub)optimal gain coefficients
%       options - struct with additional outputs, e.g. reduced-order dynamics if they are computed
%
%   Background: Computes a degree d polynomial approximation to the value function
%
%          V(x) = 1/2 ( v{2}'*(x⊗x) + ... + v{d}'*(...⊗x) )
%
%   for the polynomial control-affine system
%
%    \dot{x} = Ax + F2*(x⊗x) + F3*(x⊗x⊗x) + ...
%              + Bu + G1*(x⊗u) + G2*(x⊗x⊗u) + ...
%
%   v{2} = vec(V2) = V2(:) solves the Algebraic Riccati Equation
%
%    A'*V2 + V2*A - V2*B*R^(-1)*B'*V2 + Q = 0,
%
%   and the remaining v{i} solve linear systems arising from the
%   Hamilton-Jacobi-Bellman Partial Differential Equation. Details are in [1].
%
%   TODO:
%    - Add documentation on how to use this to return an optimal control
%    - Add documentation about the cost function that is being minimized
%
%
%   Requires the following functions from the KroneckerTools repository:
%      KroneckerSumSolver
%      kronMonomialSymmetrize
%      LyapProduct
%
%   Author: Nick Corbin, UCSD
%
%   License: MIT
%
%   Reference: [1] N. A. Corbin and B. Kramer, “The polynomial-polynomial regulator:
%              computing feedback for polynomially nonlinear systems with polynomial
%              performance indexes,” 2023.
%
%  Part of the PPR repository.
%%

% Create a vec function for readability
vec = @(X) X(:); v = cell(1,degree);

%% Process inputs
if nargin < 6
    options = struct();
end

if ~isfield(options,'skipGains'); options.skipGains = false; end
if ~isfield(options,'verbose'); options.verbose = false; end
if isfield(options,'r') && options.r ~= size(f{1}, 1); useReducedOrderModel = true; else; useReducedOrderModel = false; end

% Create pointer/shortcuts for dynamical system polynomial coefficients
if iscell(f) % polynomial drift
    A = f{1};
    lf = length(f);

    if (nargin < 5)
        degree = lf+1;
    end
else
    % Need to test and debug this option
    A = f; lf = 1; % probably still need to define degree variable
    % error("ppr: Must pass in at least quadratic dynamics")
end

if iscell(g) % polynomial input
    B = g{1};
    lg = length(g) - 1;
else
    B = g;
    lg = 0;
    g = {B};
end

[n, m] = size(B); % B is n-by-m

if iscell(q) % Polynomial state penalty
    if length(q{2}) > 1
        Q = reshape(q{2}, n, n);
    else
        Q = q{2};
    end
    lq = length(q) - 1;
else
    % quadratic state penalty with weighting matrix Q
    Q = q;
    lq = 1;
    q = {[], vec(Q)};
end

if iscell(r) % Polynomial control penalty
    if length(r{1}) > 1
        R = reshape(r{1}, m, m);
    else
        R = r{1};
    end
    lr = length(r);
else
    % quadratic control penalty with weighting matrix R
    if length(r) == 1
        r = r*eye(m);
    end
    R = r;
    lr = 1;
    r = {vec(R)};
end

% Check if R is positive definite, negative definite, or zero
if length(R) == 1
    if R > 0
        % disp('R is positive definite.')
        RPosDef = 1; % Case 1
        Rinv = 1 / R;
    elseif R < 0
        % disp('R is negative definite.')
        RPosDef = 2; % Case 2
        Rinv = 1 / R;
    else
        % disp('R is zero')
        RPosDef = 3; % Case 3
        Rinv = 0;
    end
else
    try
        chol(R);
        % disp('R is positive definite.')
        RPosDef = 1; % Case 1
    catch
        % disp('R is negative definite.')
        RPosDef = 2; % Case 2
    end
    Rinv = inv(R); % TODO: Yikes
end

%% V2, Degree 2 coefficient (k=2 case)
switch RPosDef
    case 1 % Positive definite R
        V2 = icare(A, B, Q, R);
    case 2 % Negative definite R
        if R == -1 && isscalar(Q) && Q == 0 % Computing open-loop controllability energy function; use lyap
            V2 = inv(lyap(A,(B*B.')));
        else % at least one case is when computing past energy function
            V2 = icare(A, B, Q, R, 'anti'); % alternatively -icare(-A, -B, Q, R);
        end
        if (isempty(V2) && options.verbose)
            warning('ppr: icare couldn''t find a stabilizing solution; trying the hamiltonian')
            [~, V2, ~] = hamiltonian(A, B, Q, R, true);
        end
    case 3  % Computing open-loop observability energy function; use lyap
        V2 = lyap(A.', Q);
end

if (isempty(V2))
    error('ppr: Linear system is not stabilizable')
end

%  Check the residual of the Riccati/Lyapunov equation
if (options.verbose)
    RES = A' * V2 + V2 * A - (V2 * B) * Rinv * (B' * V2) + Q;
    fprintf('  - The residual of the Riccati equation is %g\n', norm(RES, 'inf'));
    clear RES
end

%  Reshape the resulting quadratic coefficients
v{2} = vec(V2);
K1 = -R\B.'*V2;
K{1} = K1;

if useReducedOrderModel
    options.V2 = V2; options.q = q;
    options = getReducedOrderModel(f,g,options);

    % Now replace f,g,q with fr,gr,qr
    f = options.fr; g = options.gr; q = options.qr; n = options.r;
    A = f{1}; B = g{1};

    V2f = V2; K1f = K{1};
    V2 = options.Tib.'*V2*options.Tib;
    v{2} = vec(V2); K{1} = K1f*options.Tib;
    fprintf("complete. \n")

end
if ~isfield(options,'r')
    options.r = size(f{1}, 1);
end

%% v3-vd, Degree 3 and above coefficients (3<=k<=d cases)
if (degree > 2)
    % Set up the generalized Lyapunov solver (LHS coefficient matrix)
    [Acell{1:degree}] = deal((A + B * K{1}).');
    GaVb = memoize(@(a, b, v) g{a + 1}.' * sparse(reshape(v{b}, n, n^(b-1)))); % Memoize evaluation of GaVb

    % Compute the value function coefficients
    for k = 3:degree
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Form RHS vector 'b' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        b = 0; % Initialize b

        %%%%%%%%%%%%%%%%%%%%%%%%%%% Add drift components (F(x)) %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if lf > 1 % If we have polynomial dynamics
            iRange = 2:(k - 1); % Theoretical sum limits

            % iRange only needs to have at most the lf-1 last i's; for example, if there are only lf=3 and
            % there are only f{2} and f{3} to consider, iRange = [k-2, k-1] is sufficient (corresponding to
            % p=[2,3]). Otherwise f{p} doesn't exist and would require a bunch of zero entries in f{p} above lf.
            iRange = iRange(max(k - lf, 1):end); % Remove i's corresponding to non-existent f{p} entries

            for i = iRange
                p = k + 1 - i;
                b = b - LyapProduct(f{p}.', v{i}, i);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%% Add input components (G(x)) %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if R ~= 0
            for q_idx = 1:k-2                           % i + p + q = k + 1
                for p_idx = 0:lg
                    i = k+1 - p_idx - q_idx;     if i<2; break; end; if i==k; continue; end
                    b = b - i * vec(reshape(GaVb(p_idx, i, v), m, n^(k-q_idx)).' * K{q_idx});
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%% Add state penalty components (Q(x)) %%%%%%%%%%%%%%%%%%%%%%
        if k <= lq + 1
            b = b - q{k};
        end

        %%%%%%%%%%%%%%%%%%%%%% Add control penalty components (R(x)) %%%%%%%%%%%%%%%%%%%%%%
        for i = 1:lr
            for p_idx = 2:k-2 % start from 2 because the two p=1 and q=1 terms cancel with the two G terms
                q_idx = k - p_idx - i+1; %if q_idx > k-2; continue; end
                % b = b - vec(r{i}.' * kron(K{p_idx},K{q_idx})); % Naive way

                % Efficient way
                len = n^(p_idx+q_idx);
                for j = 1:size(r{i},2) % Can speed up with sparse r{i}, only iterate over nonzero columns using [~,cols,~] = find(r{i})
                    b([1:len]+(j-1)*len) = b([1:len]+(j-1)*len) ...
                        - vec(transpose(vec(K{q_idx}.' * reshape(r{i}(:,j),m,m) * K{p_idx})));
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%% Done with RHS! Now solve and symmetrize! %%%%%%%%%%%%%%%%%%%%%%
        [v{k}] = KroneckerSumSolver(Acell(1:k), b, k);
        [v{k}] = kronMonomialSymmetrize(v{k}, n, k);

        %% Now compute the gain coefficient
        % Compute last degree; GaVb for this hasn't been computed yet
        K{k-1} = zeros(m,n^(k-1));
        for q_idx = 0:lg
            j = k-q_idx;
            if j<2; break; end
            K{k-1} = K{k-1} - j/2*reshape(GaVb(q_idx,j, v), m, n^(k-1));
        end
        K{k-1} = R\K{k-1};
        % Need to add R terms

    end
end

if useReducedOrderModel
    % Convert energy functions and gain matrices back to full state dimension
    v{2} = vec(V2f);
    K{1} = K1f;
    for k = 3:degree
        v{k} = calTTv({options.TibInv}, k, k, v{k});
        if ~options.skipGains
            K{k-1} = calTTv({options.TibInv}, k-1, k-1, K{k-1}.').';
        end
    end
end


end


function [ft, gt] = linearTransformDynamics(f, g, T)
%transformDynamics Transforms the dynamics f, g by T.
%
%   Usage: [ft, gt] = transformDynamics(f, g, T)
%
%   Inputs:
%       f,g - cell arrays containing the polynomial coefficients for
%               the drift, input in the original coordinates.
%       T     - linear reduced transformation
%
%   Output:
%       ft,gt - cell arrays containing the polynomial coefficients
%                  for the transformed drift, input.
%
%   Background: Given a transformation T, compute the transformed dynamics.
%   TODO: Add more details here.
%
%   Authors: Nick Corbin, UCSD
%
% if ~iscell(h); h = {h}; end

[n, r] = size(T);
[~,m] = size(g{1});
% [p,~] = size(h{1});

ft = cell(size(f)); gt = cell(size(g)); %ht = cell(size(h));

for k = 1:length(f)
    ft{k} = zeros(n,r^k);
    for j = 1:n
        ft{k}(j,:) = ft{k}(j,:) + calTTv({T}, k, k, f{k}(j,:).').';
    end
    ft{k} = T\ft{k}; % Could rewrite to use TibInv
end

gt{1} = T\g{1};
for k = 2:length(g)
    error("Only implemented linear inputs so far!")
end

% for k = 1:length(h)
%     ht{k} = zeros(p,r^k);
%     for j = 1:p
%         ht{k}(j,:) = ht{k}(j,:) + calTTv({T}, k, k, h{k}(j,:).').';
%     end
% end


end

function options = getReducedOrderModel(f,g,options)
%getReducedOrderModel Get reduced order model
%
%   Usage:    options = getReducedOrderModel(options)
%
%   Inputs:
%       options -
%
%   Output:
%       options -
%
%   Background:
%   TODO: Add more details here.
%
%   Authors: Nick Corbin, UCSD
%

if isfield(options,'fr')
    % Already have ROM
    return;
end
if ~isfield(options,'method'); options.method = 'eigsOfV2'; end

switch options.method
    case 'eigsOfV2'
        [options.Tib, Xi] = eigs(options.V2,options.r);
        options.TibInv = options.Tib.';

        if options.verbose
            figure; semilogy(diag(Xi))
            hold on; xline(options.r)
            drawnow
        end

    case 'balancing'
        fprintf("Computing reduced order dynamics using linear balancing (r = %i, eta = %1.1f)... \n", options.r, options.eta)
        error("Need to complete this")
        % Define balancing transformation (internally balanced)
        % V2 = R*R.', W2 = L*L.'
        Lchol = lyapchol(A,B); % CHECK THIS
        Rchol = lyapchol(A.',C); % CHECK THIS

        [U, Xi, V] = svd(Lchol.'*Rchol); % from Theorem 2

        if options.verbose
            figure; semilogy(diag(Xi))
            hold on; xline(options.r)
            drawnow
        end
        % Truncate transformation to degree r
        Xi = Xi(1:options.r,1:options.r);
        V = V(:,1:options.r);
        U = U(:,1:options.r);

        % Define balancing transformation (internally balanced)
        Tib = V * diag(diag(Xi).^(-1/2));
        TibInv = diag(diag(Xi).^(-1/2)) * U.'*Lchol.';

    case 'precomputed'
        % Assume T is already given in options.Tib
        if ~isfield(options,'TibInv')
            options.TibInv = pinv(options.Tib); % yikes
        end
end

% Transform dynamics using the linear (reduced) transformation Tib
[options.fr, options.gr] = linearTransformDynamics(f, g, options.Tib);
for k = 2:length(options.q)
    if length(options.q{k}) == 1
        options.qr{k} = options.q{k};
    else
        options.qr{k} = calTTv({options.Tib}, k, k, options.q{k});
    end
end

end

