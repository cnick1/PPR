function [Z, K] = mess_care_lrradi(A, B, C, R, E)
% mess_care_lrradi Wrapper for calling mess_lrradi to solve continuous-time 
% Riccati equation using M-M.E.S.S. LR-ADI solver. 
%
%   [Z, K] = mess_care_lrradi(A, B, C, R, E) computes low-rank factors 
%   of the stabilizing solution X=ZZ' of the continuous-time algebraic
%   Riccati equation
%
%           A'X E + E'X A - E'X BR⁻¹B'X E + C'*C = 0 .
%
%   The key is that the Q matrix in the LQR problem is low-rank, Q = C'*C.
%   For this reason, we don't pass Q, we pass C, and C should be short/wide.
%   
%   This function only sets options for the solvers, and then it calls
%   functions from the M-M.E.S.S. package. Hence, the M-M.E.S.S. package
%   should be installed and added to path. 
%
%   All copyrights and credit for the M-M.E.S.S. go to the writers:
% Copyright (c) 2009-2025 Jens Saak, Martin Koehler, Peter Benner and others.
% All rights reserved.
% License: BSD 2-Clause License (see COPYING)
%

%% Usfs
opts = struct;
[oper, opts] = operatormanager(opts, 'default');

%% Options
ni = nargin;
no = nargout;

%% Process inputs
[n,m] = size(B);
if isscalar(R)
    R = R*speye(m);
end

%% Equation type
% else % Z*D*Z' case.
eqn.type = 'T';
eqn.A_ = A;
eqn.B = B;
eqn.C = C;
eqn.R = R;
eqn.haveE = true;
eqn.E_ = E;
% end

% Global.
opts.norm = 'fro';

% Shifts.
% n  = size(eqn.A_, 1);
s1 = min(n - 2, 50);
s2 = min(n - 2, 25);

opts.shifts.num_desired       = min(floor((s1 + s2) / 2) - 2, 25);
opts.shifts.history           = opts.shifts.num_desired * size(eqn.C, 1);
opts.shifts.method            = 'gen-ham-opti';
opts.shifts.naive_update_mode = false;

% RADI.
opts.radi.maxiter      = 300;
opts.radi.res_tol      = 1.0e-11;
opts.radi.rel_diff_tol = 0;
opts.radi.info         = 0;
opts.radi.trunc_tol    = 1.0e-13;

% Compute K and Z in Z*Z' format.
opts.radi.compute_sol_fac = true;
opts.radi.get_ZZt         = true;

%% Solve Equation
out = mess_lrradi(eqn, opts, oper);

if out.res(end) > opts.radi.res_tol
    mess_warn(opts, 'convergence', ...
              ['Convergence of solution only up to relative residual of %e!\n' ...
               'Check mess_lrnm and mess_lrradi for customizable solvers.'], ...
              out.res(end));
end

%% Prepare output
Z = out.Z;
K = -out.K;

end