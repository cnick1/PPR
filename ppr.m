function [v,K,options] = ppr(f, g, q, R, degree, options)
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
%                   • f must contain at least linear and quadratic coefficients
%                   • g must contain at least a linear input (B matrix)
%       q       - cell arrays containing the polynomial coefficients for the
%                 state-dependent weighting. At minimum q{2} = Q(:) where Q is
%                 the quadratic weighting, i.e. x.'*Q*x (can optionally pass Q
%                 matrix directly in this case)
%       R       - quadratic weighting matrix on the input cost
%       degree  - desired degree of the computed value function. A degree d
%                 energy function uses information from f,g,q up-to degree d-1.
%                 The default choice of d is lf+1, where lf is the degree of
%                 the drift.
%       options - struct containing additional options:
%           - verbose: print runtime information
%           - skipGains: skip computing the feedback gain coefficients
%                       (defaults to false)
%         * The next set of options are used for accelerating the computations.
%           We can use linear balanced truncation to compute a reduced-order model
%           of dimension r, and then use the reduced dynamics to compute the energy
%           functions in r dimensions rather than n; this can present a huge
%           computational speedup.
%           - r: dimension "r" of reduced order model. (Defaults to n)
%           - fr,gr,hr,qr: reduced dynamics; if for example the past energy
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
vec = @(X) X(:);

%% Process inputs
if nargin < 6
    options = struct();
end

if ~isfield(options,'skipGains'); options.skipGains = false; end
if ~isfield(options,'verbose'); options.verbose = false; end

if isfield(options,'fr')
    useReducedOrderModel = true;
    nFOM = size(f{1}, 1); % Stash original FOM n for use later
    if ~isfield(options,'r')
        options.r = size(options.fr{1}, 1);
    end
    
    % Now replace f,g,q,R with fr,gr,qr,Rr
    f = options.fr; g = options.gr;
    q = options.qr; R = options.Rr;
    
elseif isfield(options,'r') && options.r ~= size(f{1}, 1)
    useReducedOrderModel = true;
    if ~isfield(options,'eta'); options.eta = 0; end
    fprintf("Computing reduced order dynamics using linear balancing (r = %i, eta = %1.1f)... \n", options.r, options.eta)

    %% Compute reduced dynamics
    method = 2; 
    switch method
        case 1 % Use ppr recursively; requires some un-needed inversions that cause trouble sometimes (for non-minimal models)
            % First compute just the quadratic components
            v = approxPastEnergy(f, g, options.h, options.eta, 2, options.verbose);
            w = approxFutureEnergy(f, g, options.h, options.eta, 2, options.verbose);

            % Compute reduced (linear) balancing transformation and transform the dynamics
            nFOM = size(f{1}, 1); % Stash original FOM n for use later
            V2 = reshape(v{2}, nFOM, nFOM); W2 = reshape(w{2}, nFOM, nFOM);

            try
                Rchol = chol(V2, 'lower'); Lchol = chol(W2, 'lower'); % V2 = R*R.', W2 = L*L.'
            catch
                warning("ppr: Gramians not positive definite, trying sqrtm()")
                Rchol = sqrtm(V2); Lchol = sqrtm(W2); % V2 = R*R.', W2 = L*L.'
            end

            [U, Xi, V] = svd(Lchol.' / Rchol.'); % from Theorem 2
        case 2
            Pi = icare(f{1}.', options.h.', g*g.', 1/options.eta);
            P = icare(f{1}, g, options.h.'*options.h, 1/options.eta);
            try
                Rchol = chol(Pi, 'lower'); Lchol = chol(P, 'lower'); % V2 = R*R.', W2 = L*L.'
            catch
                warning("ppr: Gramians not positive definite, trying sqrtm()")
                Rchol = sqrtm(Pi); Lchol = sqrtm(P); % V2 = R*R.', W2 = L*L.'
            end
            [U, Xi, V] = svd(Lchol.' * Rchol); % from Theorem 2
    end

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
    Tib = R.' \ V * diag(diag(Xi).^(-1/2));
    TibInv = diag(diag(Xi).^(-1/2)) * U.'*Lchol.';
    
    
    % Transform dynamics using the linear (reduced) transformation Tin
    if ~iscell(g); g = {g}; end; if ~iscell(options.h); options.h = {options.h}; end
    [options.fr, options.gr, options.hr] = linearTransformDynamics(f, g, options.h, Tib);
    for k = 2:length(q)
        if length(q{k}) == 1
            options.qr{k} = q{k};
        else 
            options.qr{k} = calTTv({Tib}, k, k, q{k});
        end
    end

    clear v w V2 W2 Rchol Lchol U Xi V
    % Now replace f,g,q with fr,gr,qr
    f = options.fr; g = options.gr; q = options.qr;
    
else
    useReducedOrderModel = false;
end
if ~isfield(options,'r')
    options.r = size(f{1}, 1);
end

% Create pointer/shortcuts for dynamical system polynomial coefficients
if iscell(f)
    % polynomial drift
    A = f{1};
    % N = f{2}; % maybe don't do this here? Well if N is missing the code would break anyways
    lf = length(f);
    
    if (nargin < 5)
        degree = lf+1;
    end
else
    error("ppr: Must pass in at least quadratic dynamics")
end

if iscell(g)
    % polynomial input
    B = g{1};
    lg = length(g) - 1;
else
    B = g;
    lg = 0;
    g = {B};
end

n = size(A, 1); % A should be n-by-n
m = size(B, 2); % B should be n-by-m

if iscell(q)
    % State dependent weighting (polynomial)
    if length(q{2}) > 1
        Q = reshape(q{2}, length(A), length(A));
    else
        Q = q{2};
    end
    lq = length(q) - 1;
else
    % constant state weighting
    Q = q;
    lq = 1;
    q = {[], vec(Q)};
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
            V2 = -icare(A, B, Q, R, 'anti'); % alternatively icare(-A, -B, Q, R);
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
K{1} = -R\B.'*V2;

%% v3-vd, Degree 3 and above coefficients (3<=k<=d cases)
if (degree > 2)
    % Pre-allocate G_a.'*V_b, etc for all the a,b we need.
    % They will be computed and stored once as needed
    % For G_a, we need lg + 1 because lg is indexed from 0 (B = G0);
    % only need up to degree-1 for V_b because this will be the rhs of the last computed coefficients
    GaVb = cell(lg + 1, degree - 1);
    % GaVb = cell(2 * lg + 1, degree - 1); % For G_a, we need 2*lg + 1 why???
    
    % Compute first one
    % GaVb{1, 2} = B.' * V2;
    
    % Set up the generalized Lyapunov solver (LHS coefficient matrix)
    [Acell{1:degree}] = deal((A - (B * Rinv * B.') * V2).');
    
    % Compute the coefficients
    for k = 3:degree
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Form RHS vector 'b' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        b = 0; % Initialize b
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Add polynomial drift components %%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Add polynomial input components %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % First pre-compute GaVb components we need once to re-use without recomputing
        % Only need k-1 for V because the previous iteration for k handled k-2, etc.
        for a = 0:lg
            GaVb{a + 1, k - 1} = g{a + 1}.' * sparse(reshape(v{k - 1}, n, n ^ (k - 2)));
        end
        
        % Next add BB.' terms; this must be done separately since some of the o=0 terms need to go on the LHS
        for i = 3:k - 1 % Note how this is not 2:k; need to remove the V2 and Vk components that go in LHS
            j = k + 2 - i;
            b = b + 0.25 * i * j * vec(GaVb{1, i}.' * Rinv.' * GaVb{1, j});
        end
        
        if R ~= 0
            % Then add remaining G terms (o = 1 to 2 ell)
            for o = 1:2 * lg
                for p_idx = max(0, o - lg):min(o, lg) % These max/mins get around the nonexistent Gs, such as g{2*lg}
                    for i = 2:k - o
                        q_idx = o - p_idx;
                        j = k - o - i + 2;
                        tmp = kron(speye(n ^ p_idx), vec(speye(m)).') ...
                            * kron(vec(GaVb{q_idx + 1, j}).', kron(GaVb{p_idx + 1, i}, Rinv*speye(m))) ...
                            * kron(speye(n ^ (j - 1)), kron(perfectShuffle(n ^ (i - 1), n ^ q_idx * m), speye(m))) ...
                            * kron(speye(n ^ (k - p_idx)), vec(speye(m)));
                        b = b + 0.25 * i * j * vec(tmp);
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%% Add polynomial q(x) weighting components %%%%%%%%%%%%%%%%%%%%%%
        if k <= lq + 1
            b = b - q{k};
        end
        
        %%%%%%%%%%%%%%%%%%%%%% Done with RHS! Now solve and symmetrize! %%%%%%%%%%%%%%%%%%%%%%
        [v{k}] = KroneckerSumSolver(Acell(1:k), b, k);
        [v{k}] = kronMonomialSymmetrize(v{k}, n, k);
        
    end
    
    if ~options.skipGains
        %% Now compute the gain coefficient
        for k=2:degree-2
            K{k} = zeros(m,n^k);
            % for i=2:k % in terms of i; issues because of nonexistent Gp coefficients
            %     K{k-1} = K{k-1} - i/2*reshape(GaVb{k-i+1,i},m,n^(k-1))/R;
            % end
            % for p_idx = max(0, k-1 - lg):min(k-1, lg)
            for p_idx = 0:lg
                i = k-p_idx+1;
                if i<2; break; end
                K{k} = K{k} - transpose(i/2*reshape(GaVb{p_idx+1,i},m,n^k).'/R);
            end
        end
        
        % Compute last degree; GaVb for this hasn't been computed yet
        k = degree - 1;
        K{k} = zeros(m,n^k);
        % for i=2:k % in terms of i; issues because of nonexistent Gp coefficients
        %     K{k-1} = K{k-1} - i/2*reshape(GaVb{k-i+1,i},m,n^(k-1))/R;
        % end
        % for p_idx = max(0, k-1 - lg):min(k-1, lg)
        for p_idx = 0:lg
            i = k-p_idx+1;
            if i<2; break; end
            K{k} = K{k} - transpose(i/2*reshape(g{p_idx+1}.' * sparse(reshape(v{i}, n, n ^ (i - 1))),m,n^k).'/R);
        end
    else
        K = [];
    end
end

if useReducedOrderModel
    % Convert energy functions and gain matrices back to full state dimension
    for k = 2:degree
        v{k} = calTTv({TibInv}, k, k, v{k});
        K{k-1} = calTTv({TibInv}, k-1, k-1, K{k-1}.').';
    end
end

end


function [ft, gt, ht] = linearTransformDynamics(f, g, h, T)
%transformDynamics Transforms the dynamics f, g, h by T.
%
%   Usage: [ft, gt, ht] = transformDynamics(f, g, h, T)
%
%   Inputs:
%       f,g,h - cell arrays containing the polynomial coefficients for
%               the drift, input, and output in the original coordinates.
%       T     - linear reduced transformation
%
%   Output:
%       ft,gt,ht - cell arrays containing the polynomial coefficients
%                  for the transformed drift, input, and output.
%
%   Background: Given a transformation T, compute the transformed dynamics.
%   TODO: Add more details here.
%
%   Authors: Nick Corbin, UCSD
%
[n, r] = size(T);
[~,m] = size(g{1});
[p,~] = size(h{1});

ft = cell(size(f)); gt = cell(size(g)); ht = cell(size(h));

for k = 1:length(f)
    ft{k} = zeros(n,r^k);
    for j = 1:n
        ft{k}(j,:) = ft{k}(j,:) + calTTv({T}, k, k, f{k}(j,:).').';
    end
    ft{k} = T\ft{k};
end

gt{1} = T\g{1};
for k = 2:length(g)
    error("Only implemented linear inputs so far!")
end

for k = 1:length(h)
    ht{k} = zeros(p,r^k);
    for j = 1:p
        ht{k}(j,:) = ht{k}(j,:) + calTTv({T}, k, k, h{k}(j,:).').';
    end
end


end