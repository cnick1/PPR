function [v] = pqr(f, g, q, R, degree, verbose)
%pqr  Compute a polynomial approximation to the value function for a polynomial
% control-affine dynamical system.
%
%   Usage: v = pqr(f, g, q, R, degree, verbose)
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
%       degree  - desired degree of the computed value function. A degree d
%                 energy function uses information from f,g,q up-to degree d-1.
%                 The default choice of d is lf+1, where lf is the degree of
%                 the drift.
%       verbose - optional argument to print runtime information
%
%   Output:
%       v       - cell array containing the polynomial value function coefficients
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
%  Part of the PQR repository.
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
        degree = lf;
    end
else
    error("Must pass in at least quadratic dynamics")
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

% Check if R is positive or negative definite
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
    end
else
    try chol(R)
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
        [V2] = icare(A, B, Q, R);

        if (isempty(V2) && verbose)
            warning('pqr: icare couldn''t find stabilizing solution')
        end
    case 2 % Negative definite R
        [V2] = icare(A, B, Q, R, 'anti');

        if (isempty(V2) && verbose)
            warning('pqr: icare couldn''t find a stabilizing solution; trying the hamiltonian')
            [~, V2, ~] = hamiltonian(A, B, Q, R, true);
        end
    case 3 % Open-loop
        [V2] = lyap(A.', Q);
end

if (isempty(V2))
    error('pqr: Linear system is not stabilizable')
end

%  Check the residual of the Riccati/Lyapunov equation
if (verbose)
    RES = A' * V2 + V2 * A - (V2 * B) * Rinv * (B' * V2) + Q;
    fprintf('The residual of the Riccati equation is %g\n', norm(RES, 'inf'));
    clear RES
end

%  Reshape the resulting quadratic coefficients
v{2} = vec(V2);

%% v3, Degree 3 coefficient (k=3 case)
if (degree > 2)
    GaWb = cell(2 * lg + 1, degree - 1); % Pre-compute G_a.'*W_b, etc for all the a,b we need
    GaWb{1, 2} = B.' * V2;
    % set up the generalized Lyapunov solver
    [Acell{1:degree}] = deal((A - (B * Rinv * B.') * V2).');

    b = -LyapProduct(f{2}.', v{2}, 2);

    if lg > 0 % only for polynomial input
        Im = speye(m);
        GaWb{2, 2} = g{2}.' * V2;
        b = b + 2 * kron(speye(n ^ 3), vec(Im).') * vec(kron(GaWb{2, 2}, Rinv.' * GaWb{1, 2}));
    end

    % q(x)
    if 3 <= lq + 1
        b = b - q{3};
    end

    [v{3}] = KroneckerSumSolver(Acell(1:3), b, 3); % Solve Ax=b for k=3
    [v{3}] = kronMonomialSymmetrize(v{3}, n, 3); % Symmetrize v3

    %% k>3 cases (up to d)
    for k = 4:degree
        GaWb{1, k - 1} = B.' * reshape(v{k - 1}, n, n ^ (k - 2));

        b = -LyapProduct(f{2}.', v{k - 1}, k - 1); % Pre-compue all the L(N') terms

        % New for polynomial drift f(x)

        iRange = 2:(k - 2);
        iRange = iRange(max(k - lf, 1):end); % Need to only do lf last i's; if there are only 2 Ns for example, only need k-2! otherwise f(xi) doesnt exist and would require a bunch of empty f(xi)s
        %         TODO: may be better to write directly in terms of xi, but then how to get rid of i=k-1...?
        for i = iRange % would be from 2:k-1 but k-1 was covered in instantiation of b
            xi = k + 1 - i;
            b = b - LyapProduct(f{xi}.', v{i}, i);
        end

        % Now add all the terms from the 'B' sum by looping through the i and j
        for i = 3:(k + 1) / 2 % i+j=k+2
            j = k + 2 - i;
            tmp = GaWb{1, i}.' * Rinv.' * GaWb{1, j};
            b = b + 0.25 * i * j * (vec(tmp) + vec(tmp.'));
        end

        if ~mod(k, 2) % k is even
            i = (k + 2) / 2;
            j = i;
            tmp = GaWb{1, i}.' * Rinv.' * GaWb{1, j};
            b = b + 0.25 * i * j * vec(tmp);
        end

        % Now add the higher order polynomial terms "G" by iterating through the sums
        [g{lg + 2:2 * lg + 1}] = deal(0); % Need extra space in g because of GaWb indexing

        for o = 1:2 * lg
            for idx = 2:k - 1 % Might be repetitive
                if o + 1 < lg + 2
                    GaWb{o + 1, idx} = g{o + 1}.' * sparse(reshape(v{idx}, n, n ^ (idx - 1)));
                else
                    GaWb{o + 1, idx} = 0;
                end
            end
            for p_idx = max(0, o - lg):min(o, lg)

                for i = 2:k - o
                    q_idx = o - p_idx;
                    j = k - o - i + 2;
                    tmp = kron(speye(n ^ p_idx), vec(Im).') ...
                        * kron(vec(GaWb{q_idx + 1, j}).', kron(GaWb{p_idx + 1, i}, Rinv)) ...
                        * kron(speye(n ^ (j - 1)), kron(perfectShuffle(n ^ (i - 1), n ^ q_idx * m), Im)) ...
                        * kron(speye(n ^ (k - p_idx)), vec(Im));
                    b = b + 0.25 * i * j * vec(tmp);
                end

            end

        end

        % New for polynomial output h(x)
        % TODO: use symmetry to cut in half
        %         for p_idx = (k - lq):lq % would be 1:(k-1) but need to truncate only to h{} terms which exist
        %             q_idx = k - p_idx;
        %             b = b - vec(q{p_idx}.' * q{q_idx});
        %         end

        % q(x)
        if k <= lq + 1
            b = b - q{k};
        end

        % Done with RHS! Now solve and symmetrize!
        [v{k}] = KroneckerSumSolver(Acell(1:k), b, k);
        [v{k}] = kronMonomialSymmetrize(v{k}, n, k);
    end

end

end
