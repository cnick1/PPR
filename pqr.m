function [v] = pqr(f, g, q, R, d, verbose)
%  Calculates a polynomial approximation to the future energy function
%  for a quadratic drift, polynomial input system. The default usage is
%
%  w = pqr(f,g,h,eta,d,verbose)
%
%  where 'verbose' is an optional argument. If the system has a constant
%  input vector field Bu, the matrix B may be passes in place of a cell
%  array 'g'. The cell array 'g' should be g{1} = B = G_0, g{2} = G1, g{3} = G2 ...
%
%  Computes a degree d polynomial approximation to the future energy function
%
%          E^+(x) = 1/2 ( w{2}'*kron(x,x) + ... + w{d}'*kron(.. x) )
%
%  for the polynomial system
%
%    \dot{x} = Ax + Bu + F2*kron(x,x) + G1*kron(x,u) + G2*kron(x,x,u) + ...
%          y = Cx
%
%  where eta = 1-gamma^(-2), gamma is the parameter in the algebraic Riccati
%  equation
%
%    A'*W2 + W2*A - eta*W2*B*B'*W2 + C'*C = 0,
%
%  and in the subsequent linear systems arising from the Future H-infinity
%  Hamilton-Jacobi-Bellman Partial Differential Equation.
%
%  Note that w{2} = vec(W2) = W2(:).  Details are in Section III.B of reference [1].
%
%  Requires the following functions from the KroneckerTools repository:
%      KroneckerSumSolver
%      kronMonomialSymmetrize
%      LyapProduct
%
%  Authors: Nick Corbin, UCSD
%
%  License: MIT
%
%  Reference: [2] Scalable Computation of ℋ∞ Energy Functions for
%             Polynomial Control-Affine Systems, N. Corbin and B. Kramer,
%             arXiv:
%
%             See Algorithm 1 in [1].
%
%  Part of the NLbalancing repository.
%%

% Create a vec function for readability
vec = @(X) X(:);

%% Process inputs
if (nargin < 7)
    verbose = false;
end

% Create pointer/shortcuts for dynamical system polynomial coefficients
% polynomial drift balancing
A = f{1};
lf = length(f);

if iscell(g)
    % QB or polynomial input balancing
    B = g{1};
    lg = length(g) - 1;
else
    % Will reduce to Jeff's original code
    B = g;
    lg = 0;
    g = {B};
end

n = size(A, 1); % A should be n-by-n
m = size(B, 2); % B should be n-by-m

if iscell(q) % state dependent weighting (polynomial)
    Q = reshape(q{2},length(A),length(A));
    lq = length(q) - 1;
else % constant state weighting
    Q = q;
    lq = 1;
    q = {[],vec(Q)};
end

if length(R) == 1
    R = R * eye(m);
    Rinv = 1/R * eye(m);
else
    Rinv = inv(R);
end




%% k=2 case
% if (eta > 0) TODO: Fix various eta cases for general R
    [V2] = icare(A, B, Q, R);
    
    if (isempty(V2) && verbose)
        warning('pqr: icare couldn''t find stabilizing solution')
    end
    
% elseif (eta < 0)
%     [V2] = icare(A, B, Q, R, 'anti');
%     
%     if (isempty(V2) && verbose)
%         warning('pqr: icare couldn''t find stabilizing solution')
%         warning('pqr: using the hamiltonian')
%         [~, V2, ~] = hamiltonian(A, B, Q, R, true);
%     end
%     
% else % eta==0
%     [V2] = lyap(A.', Q);
%     
% end

if (isempty(V2))
    error('pqr: Can''t find a stabilizing solution')
end

%  Check the residual of the Riccati/Lyapunov equation
if (verbose)
    RES = A' * V2 + V2 * A - (V2 * B) * Rinv * (B' * V2) + Q;
    fprintf('The residual of the Riccati equation is %g\n', norm(RES, 'inf'));
    clear RES
end

%  Reshape the resulting quadratic coefficients
v{2} = vec(V2);

%% k=3 case
if (d > 2)
    GaWb = cell(2 * lg + 1, d - 1); % Pre-compute G_a.'*W_b, etc for all the a,b we need
    GaWb{1, 2} = B.' * V2;
    % set up the generalized Lyapunov solver
    [Acell{1:d}] = deal((A - (B * Rinv * B.')*V2).');
    
    b = -LyapProduct(f{2}.', v{2}, 2);
    
    if lg > 0 % only for polynomial input
        Im = speye(m);
        GaWb{2, 2} = g{2}.' * V2;
        b = b + 2 * kron(speye(n ^ 3), vec(Im).') * vec(kron(GaWb{2, 2}, Rinv.'*GaWb{1, 2}));
    end
    
    % q(x)
    if 3 <= lq + 1
        b = b - q{3};
    end
    
    [v{3}] = KroneckerSumSolver(Acell(1:3), b, 3); % Solve Ax=b for k=3
    [v{3}] = kronMonomialSymmetrize(v{3}, n, 3); % Symmetrize w3
    
    %% k>3 cases (up to d)
    for k = 4:d
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

% if verbose % Check HJB Residual
%     % Compute "residual":
%     %     nX = 301; nY = nX;
%     %     xPlot = linspace(-1, 1, nX);
%     %     yPlot = linspace(-1, 1, nY);
%     %     [X, Y] = meshgrid(xPlot, yPlot);
%     %     RES = zeros(nY, nX);
%     %     degree = length(w);
%     %     for i = 1:nY
%     %         for j = 1:nX
%     %             x = [X(i, j); Y(i, j)];
%     %
%     %             RES(i, j) = (0.5 * kronPolyDerivEval(w, x, degree)) * (f{1} * x + f{2} * kron(x, x)) ...
%     %                 - eta / 2 * 0.25 * kronPolyDerivEval(w, x, degree) * (g{1} + g{2} * x) * (g{1} + g{2} * x).' * kronPolyDerivEval(w, x, degree).' ...
%     %                 + 0.5 * (C * x).' * (C * x);
%     %         end
%     %     end
%     % RES = sum(sum(abs(RES))) / (nX * nY);
%     RES = computeResidualFutureHJB(f, g, h, eta, w);
%
%     fprintf('The residual of the HJB equation on the unit square is %g\n', norm(RES, 'inf'));
% else
%     RES = [];
% end

end
