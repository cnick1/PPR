function [v,K] = ppr(f, g, q, R, degree, verbose)
%ppr  Compute a polynomial approximation to the value function for a polynomial
% control-affine dynamical system.
%
%   Usage: v = ppr(f, g, q, R, degree, verbose)
% 
%       H∞ balancing energy functions can be computed as
%           [v] = ppr(f, g, cellfun(@(x) x * (-eta), h2q(h), 'un', 0), -1, degree, verbose);
%           [w] = ppr(f, g, h2q(h), 1/eta, degree, verbose);
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
%       verbose - optional argument to print runtime information
%
%   Output:
%       v       - cell array containing the polynomial value function coefficients
%       k       - cell array containing the polynomial (sub)optimal gain coefficients
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
if (nargin < 6)
    verbose = false;
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
        else
            V2 = icare(A, B, Q, R, 'anti');
        end
        if (isempty(V2) && verbose)
            warning('ppr: icare couldn''t find a stabilizing solution; trying the hamiltonian')
            [~, V2, ~] = hamiltonian(A, B, Q, R, true);
        end
    case 3 % Open-loop
        V2 = lyap(A.', Q);
end

if (isempty(V2))
    error('ppr: Linear system is not stabilizable')
end

%  Check the residual of the Riccati/Lyapunov equation
if (verbose)
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

    %% Now compute the gain coefficient
    for k=2:degree-2
        K{k} = zeros(m,n^k); 
        % for i=2:k % in terms of i; issues because of nonexistent Gp coefficients
        %     K{k-1} = K{k-1} - i/2*reshape(GaVb{k-i+1,i},m,n^(k-1))/R;
        % end
        % for p_idx = max(0, k-1 - lg):min(k-1, lg)
        for p_idx = 0:lg
            i = k-p_idx+1; 
            if i<2; break; end;
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
        if i<2; break; end;
        K{k} = K{k} - transpose(i/2*reshape(g{p_idx+1}.' * sparse(reshape(v{i}, n, n ^ (i - 1))),m,n^k).'/R);
    end
end

end
