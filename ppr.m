function [v,K,options] = ppr(f, g, q, r, degree, options)
%ppr Compute a polynomial approximation to the optimal feedback law
%    and value function for a polynomial control-affine dynamical system.
%
%   Usage: [v, K] = ppr(f, g, q, r, degree)
%
%       PPR controller can be computed as [~, K] = ppr(f, g, q, R, degree)
%           Recommended implementation of PPR controllers for simulations
%           is to define the dynamics and controller using the kronPolyEval
%           function:
%               uPPR = @(x) kronPolyEval(K, x);
%               F = @(x) kronPolyEval(f, x);
%               G = @(x) kronPolyEval(g, x, scenario='G(x)');
%           The system can then be simulated for example using ode45:
%               [t, X] = ode45(@(t, x) F(x) + G(x) * uPPR(x), tspan, x0);
%           In some cases, it is more efficient or more accurate to program
%           F(x) and G(x) differently, e.g. instead of using the polynomial
%           approximation given by the cell arrays f,g, sometimes the full
%           nonlinear functions are known, e.g. sin(x), etc.
%
%       H∞ balancing energy functions can be computed as
%           [v] = ppr(f, g, cellfun(@(x) x * (-eta), h2q(h), 'un', 0), -1, degree);
%           [w] = ppr(f, g, h2q(h), 1/eta, degree);
%
%       For dynamics in descriptor form with an invertible mass matrix E,
%           sparsity can be leveraged by passing the mass matrix E as an
%           optional argument.
%               options.E = [mass matrix];
%               [v, K] = ppr(f, g, q, r, degree, options)
%           In this case, the value function is computed for a transformed
%           coordinate system z = E*x, so to evaluate it one would call
%               V = @(x) kronPolyEval(v, E*x);
%
%   Inputs:
%       f,g     - cell arrays containing the polynomial coefficients
%                 for the drift and input.
%       q       - cell arrays containing the polynomial coefficients for the
%                 state penalty in the cost. At minimum q{2} = Q(:) where Q is
%                 the quadratic weighting, i.e. xᵀ Q x (can optionally pass Q
%                 matrix directly in this case)
%       r       - cell arrays containing the polynomial coefficients for the
%                 control penalty in the cost. At minimum r{1} = R where R is
%                 the quadratic weighting, i.e. uᵀ R u (can optionally pass R
%                 matrix directly in this case)
%       degree  - desired degree of the computed value function. A degree d
%                 energy function uses information from f,g,q up-to degree d-1.
%                 The default choice of d is lf+1, where lf is the degree of
%                 the drift.
%       options - struct containing additional options:
%           • E: nonsingular mass matrix, for dynamics Eẋ = f(x) + g(x)u
%           • verbose: print runtime information
%         * The next set of options are used for accelerating the computations.
%           We can use linear balanced truncation to compute a reduced-order model
%           of dimension r, and then use the reduced dynamics to compute the energy
%           functions in r dimensions rather than n; this can present a huge
%           computational speedup.
%           • r: dimension "r" of reduced order model. (Defaults to n)
%           • fr,gr,qr: reduced dynamics; if for example the past energy
%             function has already been computed, instead of recomputing the
%             reduced dynamics, they can be passed in directly.
%           • eta: H-infinity balancing parameter. Defaults to open-loop balancing
%                 (eta = 0).
%
%   Output:
%       v       - cell array containing the polynomial value function coefficients
%       K       - cell array containing the polynomial optimal gain coefficients
%       options - struct with additional outputs, e.g. reduced-order dynamics if they are computed
%
%
%   Description: We seek a solution to the optimal control problem
%
%           minᵤ   J(x,u) = ½∫ xᵀ Q(x) x + uᵀ R(x) u dt
%           s.t.    ẋ = f(x) + g(x) u,   x(0) = x₀
%
%   The solution is given by the solution to the HJB PDEs
%
%       (1)  0 = 𝜕ᵀV(x)/𝜕x (f(x) + g(x) u) + ½ xᵀ Q(x) x + ½ uᵀ(x) R(x) u(x)
%       (2)  0 = R(x) u(x) + gᵀ(x) 𝜕V(x)/𝜕x
%
%   In theory, the second HJB PDE can be solved for the optimal control in
%   terms of the gradient of the value function, which can then be used to
%   eliminate u(x) from the first HJB PDE so that it becomes independent of
%   u(x) and only depends on V(x):
%
%        u(x) = -R⁻¹(x) gᵀ(x) 𝜕V(x)/𝜕x
%           0 = 𝜕ᵀV(x)/𝜕x f(x) - ½ 𝜕ᵀV(x)/𝜕x g(x) R⁻¹(x) gᵀ(x) 𝜕V(x)/𝜕x + ½ xᵀ Q(x) x
%
%   However, it is advantageous numerically to treat the two equations
%   separately, alternating successively solving a term in one and then the
%   other. We compute solution approximations to the HJB PDEs (1) and (2)
%   using the method of Al'brekht [2], i.e. we compute the Taylor expansions:
%
%           V(x) = 1/2 ( v₂ᵀ(x ⊗ x) + v₃ᵀ(x ⊗ x ⊗ x) + ... +   vᵈᵀ(... ⊗ x) )
%           u(x) = K₁ x + K₂(x ⊗ x) +  K₃(x ⊗ x ⊗ x) + ... +  Kᵈ⁻¹(... ⊗ x) )
%
%   based on the Taylor expansions for the dynamics, written as
%
%           ẋ = A x + F₂ (x ⊗ x) + F₃ (x ⊗ x ⊗ x) + ...
%               + B u + G₁ (x ⊗ u) + G₂ (x ⊗ x ⊗ u) + ...
%
%   and the Taylor expansion of the cost, written in a similar form as
%
%           J(x,u) = ½∫ xᵀ Q x + uᵀ R u + q₃ᵀ(x ⊗ x ⊗ x) + ...
%                       + xᵀr₁ᵀ(u ⊗ u) + ... dt
%
%   Inserting all these polynomial expressions into the HJB PDEs (1) and (2)
%   leads to equations for the value function coefficients v₂, v₃,..., vᵈ
%   and the optimal feedback law gain coefficients K₁, K₂, K₃,..., Kᵈ⁻¹.
%   The first coefficients are the LQR solutions; v₂ = vec(V₂) = V₂(:)
%   solves the Algebraic Riccati Equation
%
%           Aᵀ V₂ + V₂ A - V₂ B R⁻¹ Bᵀ V₂ + Q = 0,
%
%   and K₁ = -R⁻¹ Bᵀ V₂. The remaining vᵢ solve linear systems arising
%   from (1), and the remaining Kᵢ solve linear systems arising from
%   (2). The details can be found in [1]. The theoretical results can be
%   found in [2,3]. This work is heavily inspired by and builds on the
%   works [4,5,6].
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
%   Reference: [1] N. A. Corbin and B. Kramer, "Computing solutions to the
%               polynomial-polynomial regulator problem,” in 2024 63rd IEEE
%               Conference on Decision and Control, Dec. 2024
%              [2] E. G. Al’brekht, “On the optimal stabilization of
%               nonlinear systems," Journal of Applied Mathematics and
%               Mechanics, vol. 25, no. 5, pp. 1254–1266, Jan. 1961, doi:
%               10.1016/0021-8928(61)90005-3
%              [3] D. L. Lukes, “Optimal regulation of nonlinear dynamical
%               systems,” SIAM Journal on Control, vol. 7, no. 1, pp.
%               75–100, Feb. 1969, doi: 10.1137/0307007
%              [4] A. Krener, C. Aguilar, and T. Hunt, “Series solutions of
%               HJB equations,” Mathematical System Theory—Festschrift in
%               Honor of Uwe Helmke on the Occasion of his Sixtieth
%               Birthday, pp. 247–260, 2013
%              [5] J. Borggaard and L. Zietsman, “The quadratic-quadratic
%               regulator problem: approximating feedback controls for
%               quadratic-in-state nonlinear systems,” in 2020 American
%               Control Conference (ACC), Jul. 2020, pp. 818–823. doi:
%               10.23919/ACC45564.2020.9147286
%              [6] J. Borggaard and L. Zietsman, “On approximating
%               polynomial-quadratic regulator problems,” IFAC-PapersOnLine,
%               vol. 54, no. 9, pp. 329–334, 2021, doi:
%               10.1016/j.ifacol.2021.06.090.
%
%  Part of the PPR repository.
%
%  See also: KroneckerSumSolver, kronMonomialSymmetrize, LyapProduct
%%

% Create a vec function for readability
vec = @(X) X(:); v = cell(1,degree); K = cell(1,degree-1);

%% Process inputs
if nargin < 6
    options = struct();
end

if ~isfield(options,'verbose'); options.verbose = false; end
if isfield(options,'reducedDimension') && options.reducedDimension ~= size(f{1}, 1); useReducedOrderModel = true; else; useReducedOrderModel = false; end
if ~isfield(options,'E'); options.E = []; end
if ~isfield(options,'solver'); options.solver = []; end
E = options.E;  % full?

% Create pointer/shortcuts for dynamical system polynomial coefficients
if iscell(f) % polynomial drift
    A = f{1}; % full?
    lf = length(f);
    
    if (nargin < 5)
        degree = lf+1;
    end
else
    A = f; lf = 1; % full?
    f = {A};
end

if iscell(g) % polynomial input
    B = g{1}; % full?
    lg = length(g) - 1;
else
    B = g;  % full?
    lg = 0;
    g = {B};
end

[n, m] = size(B); % B is n-by-m

if iscell(q) % Polynomial state penalty
    if length(q{2}) > 1
        Q = reshape(full(q{2}), n, n);
    else
        Q = q{2};
    end
    lq = length(q) - 1;
    
    % If q's are scalars, use 'identity' tensors
    for i=3:length(q)
        if isscalar(q{i})
            q{i} = q{i}*sparse(linspace(1,n^i,n),1,1); % Like an identity matrix but a tensor, vectorized
        end
    end
else
    % quadratic state penalty with weighting matrix Q
    Q = q;
    lq = 1;
    q = {[], vec(Q)};
end

% Check if R is positive definite, negative definite, or zero
if ~iscell(r) && isscalar(r)
    if r > 0
        % disp('R is positive definite.')
        RPosDef = 1; % Case 1
        Rinv = 1 / r;
    elseif r < 0
        % disp('R is negative definite.')
        RPosDef = 2; % Case 2
        Rinv = 1 / r;
    else
        % disp('R is zero')
        RPosDef = 3; % Case 3
        Rinv = 0;
    end
else
    try
        chol(r{1});
        % disp('R is positive definite.')
        RPosDef = 1; % Case 1
    catch
        % disp('R is negative definite.')
        RPosDef = 2; % Case 2
    end
    Rinv = inv(r{1}); % TODO: Yikes
end

if iscell(r) % Polynomial control penalty
    if length(r{1}) > 1
        R = reshape(r{1}, m, m);
    else
        R = r{1};
    end
    lr = length(r) - 1;
    
    % If r's are scalars, use 'identity' tensors
    for i=2:length(r)
        if isscalar(r{i})
            % Naive way
            % temp = sparse(m^2,n^(i-1)); temp(:,linspace(1,n^(i-1),n)) = repmat(vec(speye(m)),1,n);
            % r{i} = r{i}*temp; % Like an identity matrix but a tensor, vectorized
            
            % More efficient way: rows=linspace(1,m^2,m), columns=linspace(1,n^(i-1),n), so use
            % repmat to duplicate rows n times and columns m times, forming all possible combinations.
            r{i} = r{i}*sparse(repmat(linspace(1,m^2,m),1,n),repmat(linspace(1,n^(i-1),n),m,1),1);
        end
    end
else
    % quadratic control penalty with weighting matrix R
    if isscalar(r)
        r = r*eye(m);
    end
    R = r;
    lr = 0;
    r = {vec(R)};
end


%% V2, Degree 2 coefficient (k=2 case)
if options.verbose; T1=tic; end
switch RPosDef
    case 1 % Positive definite R
        if isfield(options,'lrradi') && options.lrradi
            try
                try
                    addpath('../mmess/'); mess_path;
                catch ME
                    % if strcmp(ME.identifier,'MESS:path_exists')
                    %     warning(ME.message)
                    % end
                end
                [Z,K1] = mess_care_lrradi(A, B, options.C, R, E);
                V2 = factoredMatrix(Z);
            catch ME
                if ~isfield(options,'C')
                    error("M-M.E.S.S. solver requires low-rank Q = C'*C; pass in C using options.C")
                else
                    error('Something went wrong calling M-M.E.S.S.; did you clone it into the parent directory?')
                end
            end
        else
            V2 = icare(full(A), full(B), full(Q), full(R), [], full(E));
        end
    case 2 % Negative definite R
        if isscalar(Q) && Q == 0 && norm(R + eye(length(R))) < 1e-12 % Computing open-loop controllability energy function; use lyap
            V2 = factoredMatrixInverse(lyapchol(A, B, E));
        else % At least one case is when computing past energy function
            V2 = icare(A, B, Q, R, [], E, 'anti'); % alternatively -icare(-A, -B, Q, R);
        end
        if (isempty(V2) && options.verbose)
            warning('ppr: icare couldn''t find a stabilizing solution; trying the hamiltonian')
            [~, V2, ~] = hamiltonian(A, B, Q, R, true);
        end
    case 3  % Computing open-loop observability energy function; use lyap
        if isa(q{2}, 'factoredMatrix')
            V2 = factoredMatrix(lyapchol(A', q{2}.Z, E').');
        else
            V2 = lyap(A', Q, [], E');
        end
end
if options.verbose; fprintf("Riccati solution computed in %2.2f seconds ... ", toc(T1)); end

if isempty(V2)
    error('ppr: Linear system is not stabilizable')
end

%  Check the residual of the Riccati/Lyapunov equation
if options.verbose
    if isempty(E); RES = A' * V2     +      V2 * A - (    V2 * B) * Rinv * (B' * V2    ) + Q;
    else;          RES = A' * (V2 * E) + E' * (V2 * A) - (E' *(V2 * B)) * Rinv * ((B' * V2) * E) + Q; end
    fprintf('  - The residual of the Riccati equation is %g\n', norm(RES, 'inf')); clear RES
end

%  Reshape the resulting quadratic coefficients
v{2} = vec(V2); % possibly factoredMatrix
if ~exist('K1','var')
    if isempty(E); K1 = -Rinv*B.'*V2;
    else;          K1 = -Rinv*B.'*V2*E; end
end
K{1} = K1;

if useReducedOrderModel
    options.V2 = V2; options.q = q; options.R = r;
    options = getReducedOrderModel(f,g,options);
    
    % Now replace f,g,q with fr,gr,qr
    f = options.fr; g = options.gr; q = options.qr; r = options.Rr; E = options.Er; n = options.reducedDimension;
    A = f{1}; B = g{1};
    
    V2f = options.V2; K1f = K{1}; % At this stage, V2 is replaced with V2bar in the E coordinates... this may be an issue down the road elsewhere
    V2 = options.T.'*(V2f*options.T); % no longer factoredMatrix
    v{2} = vec(V2); K{1} = K1f*options.T;
end
if ~isfield(options,'reducedDimension')
    options.reducedDimension = size(f{1}, 1);
end

%% v3-vd, Degree 3 and above coefficients (3<=k<=d cases)
if (degree > 2)
    % Set up the generalized Lyapunov solver (LHS coefficient matrix)
    [Acell{1:degree}] = deal((A + B * K{1}).');
    % GaVb = memoize(@(a, b, v) g{a + 1}.' * sparse(reshape(v{b}, n, n^(b-1)))); % Memoize evaluation of GaVb (I had it cast to sparse for some reason, not sure why)
    if isempty(E)
        GaVb = memoize(@(a, b, v) g{a + 1}.' * reshape(v{b}, n, n^(b-1))); % Memoize evaluation of GaVb
    else
        GaVb = memoize(@(a, b, v) g{a + 1}.' * kroneckerRight(reshape(v{b}, n, n^(b-1)),E)); % Memoize evaluation of GaVb with E matrix
    end
    
    for k = 3:degree
        %% Compute the value function coefficient
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Form RHS vector 'b' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        b = zeros(n^k,1); % Initialize b
        K{k-1} = zeros(m,n^(k-1)); % Initialize K
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Add drift components (F(x)) %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if lf > 1 % If we have polynomial drift dynamics
            % Here we add the f(x) terms. The range for possible i's is
            % 2:k-1, since vₖ is in the left-hand-side terms. However, not all of
            % these i's have corresponding Fₚ's that exist. For example, maybe only
            % F₂ and F₃ exist, in which case only the last two i's, [k-2, k-1], are
            % required. Indexing on the i's backwards and checking if we have run out
            % of Fₚ's allows us to do this easily.
            
            for i = flip(2:k-1)                           % Theoretical sum limits for Vᵢ's; index backwards in i so that p indexes forwards
                p = k + 1 - i;
                if p > lf; break; end                     % Only run while we have F_p's left
                
                b = b - LyapProduct(f{p}.', v{i}, i, E);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Add input components (G(x)) %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if RPosDef ~= 3 % as long as R is not 0
            % Here we add the g(x)u(x) terms. The sum is over indexes satisfying
            % i+p+q=k+1. The range for possible q's is 1:k-2 (technically q could go
            % up to k-1, but that term cancels with some of the R(x) terms), the range
            % for i's is 2:k-1, and the range for possible p's is 0:lg. However, the q=1
            % and p=0 term belongs on the left-hand-side, so we need to skip that term
            % with a continue statement. And since q and p index forwards, i indexes
            % backwards (at most from k), but once i reaches i<2 we run out of Vᵢ's,
            % so we need a break statement.
            
            for q_idx = 1:k-2                             % Theoretical sum limits for K_q's
                for p_idx = 0:lg                          % Theoretical sum limits for G_p's
                    i = k+1 - p_idx - q_idx;
                    if i==k; continue; end                % Skip the LHS term that is already accounted for
                    if i<2; break; end                    % Only run while we have Vᵢ's left
                    % b = b - i * vec(K{q_idx}.' * reshape(GaVb(p_idx, i, v), m, n^(k-q_idx)));
                    b = b - i * vec(K{q_idx}.' * reshape(GaVb(p_idx, i, v), m, []));
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%% Add state penalty components (Q(x)) %%%%%%%%%%%%%%%%%%%%%%
        if k <= lq + 1                                    % Simply check if q{k} exists
            b = b - q{k};
        end
        
        %%%%%%%%%%%%%%%%%%%%%% Add control penalty components (R(x)) %%%%%%%%%%%%%%%%%%%%%%
        % Here we add the u(x)ᵀR(x)u(x) terms. The sum is over indices satisfying
        % i+p+q=k. Note, the r cell array does not use zero indexing. The range for
        % possible i's is 0:lr. The theoretical range for p's and q's would be
        % 1:k-1, but the k-1 term hasn't been computed yet! It turns out these terms
        % (i=0,p=1,q=k-1 & i=0,p=k-1,q=1) cancel with terms in the G(x) sums. So we
        % will index from 1:k-2 and skip the k-1 term. (Note that p=1 and q=1 DO
        % appear if lr>0, so we can't just index 2:k-2.)
        for i = 0:lr                                     % Theoretical sum limits for r_i's
            for p_idx = 1:k-2                            % Theoretical sum limits for K_p's
                q_idx = k - p_idx - i;
                if q_idx==k-1; continue; end             % Skip the k-1 term that cancels with the G(x) terms
                if q_idx<1; break; end                   % Only run while we have K_q's left
                
                len = n^(p_idx+q_idx);
                for j = 1:size(r{i+1},2) % Can speed up with sparse r{i}, only iterate over nonzero columns using [~,cols,~] = find(r{i})
                    bRange = (1:len)+(j-1)*len;
                    b(bRange) = b(bRange) ...
                        - vec(transpose(vec(K{q_idx}.' * reshape(r{i+1}(:,j),m,m) * K{p_idx})));
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%% Done with RHS! Now symmetrize and solve! %%%%%%%%%%%%%%%%%%%%%%
        b = kronMonomialSymmetrize(b, n, k);
        [v{k}] = KroneckerSumSolver(Acell(1:k), b, k, E, [], options.solver);
        
        %% Now compute the gain coefficient
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Add input components (G(x)) %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Here we add the terms due to the gradient of the value function. The
        % indices must satisfy j+q=k. The range for the q's is 0:lg; the
        % corresponding range for the js is k:2, assuming the G's exist. In
        % other words, we start with the term given by B and the last computed
        % vk, and go down in the vks from there. We will go as high as there are
        % q's, so index forward on the q's. If you run out of V coeffients, you
        % are done and can break, even if more G's remain.
        for q_idx = 0:lg                                     % Theoretical sum limits for r_i's
            j = k - q_idx;
            if j<2; break; end                               % Only run while we have V_j's left
            K{k-1} = K{k-1} - j/2*reshape(GaVb(q_idx,j, v), m, n^(k-1));
        end
        
        %%%%%%%%%%%%%%%%%%%%% Add control penalty components (R(x)) %%%%%%%%%%%%%%%%%%%%%%%%
        % Here we add the terms due to the nonlinear control penalty R(x). The
        % theoretical range for q would be 0:lr, but the q=0 term goes on the
        % left-hand-side, so the range for q is 1:lr. Then the indices must satisfy
        % p+q=k-1.
        for q_idx = 1:lr                                     % Theoretical sum limits for r_i's
            p = k-1 - q_idx;
            if p<1; break; end                               % Only run while we have K_p's left
            Rq = reshape(r{q_idx+1}.',m*n^q_idx,m); % Reshape doesn't copy, just creates pointer
            for j=1:m
                K{k-1}(j,:) = K{k-1}(j,:) - vec(reshape(Rq(:,j),n^q_idx,m)*K{p}).';
            end
        end
        
        % Now multiply by R₀⁻¹
        K{k-1} = Rinv*K{k-1};
    end
end

if useReducedOrderModel
    % Use special data structure to avoid storing huge coefficients
    v{2} = vec(V2f);
    v = factoredValueArray(v, options.TInv);
    K{1} = K1f;
    K = factoredGainArray(K, options.TInv);
end


end

function options = getReducedOrderModel(f,g,options)
%getReducedOrderModel Apply linear model reduction to get a ROM
%
%   Usage:    options = getReducedOrderModel(f,g,options)
%
%   Inputs:
%       f,g     - cell arrays containing the polynomial coefficients
%                 for the drift and input.
%       options - struct containing additional quantities to transform:
%        • q: cell arrays containing the polynomial coefficients for the
%                 state penalty in the cost.
%        • r: cell arrays containing the polynomial coefficients for the
%                 control penalty in the cost.
%        • E: nonsingular mass matrix, for dynamics Eẋ = f(x) + g(x)u.
%
%   Output:
%       options - struct containing transformed (reduced) quantities:
%        • fr,gr: reduced drift and input.
%        • qr: reduced state penalty weights.
%        • Rr: reduced control penalty weights.
%        • Er: in the reduced form, Er=[] always.
%
%   Background: Given a transformation x = Txᵣ, we seek to represent the
%    dynamics for the control-affine system
%        Eẋ = f(x) + g(x) u
%    in the new coordinates as
%        ẋᵣ = f(xᵣ) + g(xᵣ) u
%    Applying the transformation yields
%        ETẋᵣ = f(Txᵣ) + g(Txᵣ) u
%    Here we have a choice: do we multiply by E⁻¹ and THEN enforce
%    orthogonality by Tᵀ, or do we enforce orthogonality by Tᵀ and then
%    multiply by Eᵣ⁻¹. These are not in general the same, due to the fact
%    that TᵀT = I but TTᵀ ≠ I. In this function, we do the latter option,
%    although this seems to cause issues. Multiplying by Tᵀ, we get
%        TᵀETẋᵣ = Tᵀf(Txᵣ) + Tᵀg(Txᵣ) u
%    If E=I, then this is already in standard form. Otherwise, multiplying
%    by the inverse of Eᵣ=TᵀET then puts the model in standard form:
%        ẋᵣ = fᵣ(xᵣ) + gᵣ(xᵣ) u
%    The coefficients of fᵣ(xᵣ) and gᵣ(xᵣ) in polynomial form can be
%    obtained from the coefficients of f(x) and g(x), and we have to
%    transform q(x) and r(x) similarly. The transformation is mainly
%    applied by repeated Kronecker multiplication with T:
%        Eẋ  = A x + F₂ (x ⊗ x) + ...
%       ETẋᵣ = A Txᵣ + F₂ (Txᵣ ⊗ Txᵣ) + ...
%     TᵀETẋᵣ = TᵀAT xᵣ + TᵀF₂(T⊗T) (xᵣ ⊗ xᵣ) + ...
%            = Aᵣ xᵣ + F₂ᵣ (xᵣ ⊗ xᵣ) + ...
%
%   Authors: Nick Corbin, UCSD
%
vec = @(X) X(:);

if isfield(options,'fr')
    return; % Already have ROM
end
if ~isfield(options,'method'); options.method = 'eigsOfV2'; end
if true || options.verbose; T1 = tic; end

switch options.method
    case 'eigsOfV2'
        if isa(options.V2,'factoredMatrix')
            if ~isempty(options.E)
                options.V2.Z = options.E.'*options.V2.Z;
            end
            [options.T, Xi, ~] = svd(options.V2.Z, 'econ');
            Xi = Xi(1:options.reducedDimension,1:options.reducedDimension).^2;
            options.T = options.T(:,1:options.reducedDimension);
        else
            if ~isempty(options.E)
                options.V2 = options.E.'*options.V2*options.E;
            end
            [options.T, Xi] = eigs(full(options.V2),options.reducedDimension);
        end
        
        options.TInv = options.T.';
        if options.verbose
            figure; semilogy(diag(Xi)); hold on; xline(options.reducedDimension); drawnow
        end
        
        % Could have other cases like Balanced Truncation, POD, etc.
end

% Transform dynamics using the linear (reduced) transformation T
[n, r] = size(options.T);
[~, m] = size(g{1});

%% Transform f(x)
options.fr = cell(size(f));
if issparse(f{end})
    % Use sparse indexing to form fr more efficiently
    options.fr{1} = options.TInv*(f{1}*options.T);
    for k = 2:length(f)
        if nnz(f{k}) % skip all zero coefficients
            % Option 4: Loop over r^k columns
            TinvF = options.TInv*f{k}; % do this first to reduce the row dimension (column sparsity is still preserved)
            [Fi, Fj, Fv] = find(TinvF);
            options.fr{k} = zeros(r,r^k); % result is dense, so better to use zeros
            rowinds = cell(1, k); % Preallocate cell array for k row indices
            [rowinds{:}] = ind2sub(repmat(n, 1, k), Fj);
            
            for q = 1:r^k % due to inversion, columns are dense; could maybe do this block-wise but this works well enough since r is small (I did try but then memory becomes an issue)
                colinds = cell(1, k); % Preallocate cell array for k column indices
                [colinds{:}] = ind2sub(repmat(r, 1, k), q);
                
                % Efficient sparse evaluation qth column of f{k}*(T⊗...⊗T)
                Tprod = ones(size(Fj));
                for p = 1:k
                    Tprod = Tprod .* options.T(uint32(rowinds{p}),uint32(colinds{p})); %uint necessary because doubles are sometimes not exact integers
                end
                
                options.fr{k}(:,q) = accumarray(Fi, Fv .* Tprod, [r 1]);
            end
        else
            options.fr{k} = sparseIJV(r,r^k); % sparseIJV is best for empty arrays
        end
    end
    
else
    for k = 1:length(f)
        % options.fr{k} = options.TInv*kroneckerRight(f{k},options.T);
        options.fr{k} = kroneckerRight(options.TInv*f{k},options.T); % More efficient to multiply by Tinv first
    end
end
%% Transform g(x)
% Instead of dealing with g(x) matrix directly, deal with g_i(x) vectors so
% I can use the same code as for f(x). So loop over m input channels, use gttemp,
% and convert to gr after
lg = length(g);
gttemp = cell(lg,m);
for k = 1:m
    for i = 1:(lg-1) % this internal code block is identical to above for f(x)
        gttemp{i+1,k} = kroneckerRight(g{i+1}(:,k:m:end),options.T); % k:m:end needed for converting from g(x)u to sum g_i(x) u_i
    end
end

% Now convert from g_i(x) vectors back to g(x) matrix
options.gr{1} = options.TInv*g{1}; % B is not state dependent
for i=1:(lg-1)
    options.gr{i+1} = zeros(n,m*n^i);
    for kk=1:m
        options.gr{i+1}(:,kk:m:end) = gttemp{i+1,kk};
    end
    options.gr{i+1} = options.TInv*options.gr{i+1};
end

%% Transform Q(x)
for k = 2:length(options.q)
    if isscalar(options.q{k})
        options.qr{k} = options.q{k};
    else
        % options.qr{k} = kroneckerRight(options.q{k}.',options.T).';
        options.qr{k} = kroneckerLeft(options.T.',options.q{k});
    end
end

%% Transform R(x)
options.Rr{1} = options.R{1};
for k = 2:length(options.R)
    if isscalar(options.R{k})
        options.Rr{k} = options.R{k};
    else
        options.Rr{k} = kroneckerRight(options.R{k},options.T);
    end
end

%% Transform E
if ~isempty(options.E)
    %% Transform E
    options.Er = options.TInv*options.E*options.T;
    
    %% Put in standard state-space form
    % Since the reduced system is dense anyways, we can invert Er and put
    % in standard state-space form to use more efficient Lyapunov solvers
    
    % Transform f(x)
    for k = 1:length(f)
        options.fr{k} = options.Er\options.fr{k};
    end
    
    % Transform g(x)
    for k = 1:length(g)
        options.gr{k} = options.Er\options.gr{k};
    end
end
options.Er = []; % ROM is always in standard form
if true || options.verbose; fprintf("ROM computed in %2.2f seconds ... \n", toc(T1)); end
end
