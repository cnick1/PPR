function v = zubov(f, q, degree, verbose)
%zubov Estimate the region of attraction for a nonlinear system.
%
%   Usage: v = zubov(f, q, degree)
%
%   Inputs:
%       f       - cell arrays containing the polynomial coefficients
%                 for the drift
%       degree  - desired degree of the computed Lyapunov function candidate.
%                 A degree d approximation uses information from f up to
%                 degree d-1. The default choice of d is lf+1, where lf is
%                 the degree of the drift.
%       verbose - optional argument to print runtime information
%
%   Output:
%       v       - cell array containing the candidate Lyapunov function
%                 coefficients
%
%   Description: We seek to characterize the region of attraction Rₐ of the
%   equilibrium at the origin for the uncontrolled nonlinear system
%
%       (1)  ẋ = f(x)
%
%   The main theory on this topic is due to Zubov [1], building on the
%   prior related foundational work of Lyapunov [2]. A power-series
%   approach to solving the Zubov equation was detailed in [3], and the
%   topic was also discussed in [4]. My present understanding of the topic
%   goes as follows:
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %%%%%%%%%%%%%%%%%%% Approach 1: Lyapunov's Theorem %%%%%%%%%%%%%%%%%%%%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [Thm 4.1, 5] Let x = 0 be an equilibrium point of (1) and D ⊂ Rⁿ be a
%   domain containing x = 0. Let V:D → R be a continuously differentiable
%   function such that
%        (a) V(0) = 0  and  V(x) > 0  in  D − {0}
%        (b) ̇V(x) < 0  in  D − {0}
%   Then, x=0 is asymptotically stable.
%
%   [Lemma 8.1, 5] If x = 0 is an asymptotically stable equilibrium point
%   for (1), then its region of attraction Rₐ is an open, connected,
%   invariant set. Moreover, the boundary of Rₐ is formed by trajectories.
%
%   As Khalil points out, D is not Rₐ, because while D guarantees that
%   trajectories will move down the Lyapunov function, it is not
%   necessarily positively invariant, so trajectories may leave D and then
%   may not move down the Lyapunov function. What is needed is a compact
%   positively invariant subset of D; one way to construct this is with a
%   Lyapunov sublevel set:
%
%           Ωc = {x ∈ Rⁿ | V(x) ≤ c}
%
%   Although not all Ωc are subsets of D, the largest Ωc ⊂ D is an estimate
%   of the region of attraction. (Proving Ωc ⊂ D can be done only in
%   special cases, otherwise numerical evaluation may be necessary which of
%   course is not rigorous). This is positively invariant because if you
%   start in Ωc AND Ωc ⊂ D, then since ̇V(x) < 0 in D, V(x) ≤ c will always
%   hold, meaning that you stay in Ωc i.e. it is positively invariant.
%
%   So how do we construct V(x)? By solving ̇V(x) < 0. We can rewrite this
%   as 𝜕ᵀV(x)/𝜕x f(x) < 0, which can also be thought of as
%
%       (2) 𝜕ᵀV(x)/𝜕x f(x) + q(x) = 0 , q(x) > 0    e.g.   q(x) = ½ xᵀ Q(x) x
%
%   This is a nonlinear Lyapunov-type equation; the keen-eyed will notice
%   that this is in fact very similar to the observability energy HJB equation.
%   Indeed, for linear systems, we can write this as
%
%       xᵀ V₂ A x + ½ xᵀ Q x = 0 ,             Q > 0
%       xᵀ Aᵀ V₂x + xᵀ V₂ A x + xᵀ Q x = 0 ,   Q > 0
%       Aᵀ V₂ + V₂ A + Q = 0 ,                 Q > 0
%
%   For nonlinear systems, we can solve 𝜕ᵀV(x)/𝜕x f(x) + q(x) = 0 via
%   Taylor expansions; any positive definite function q(x) will work.
%   Furthermore, any approximation V(x) will work, as long as you can
%   guarantee
%        (a) V(0) = 0  and  V(x) > 0  in  D − {0}
%        (b) ̇V(x) < 0  in  D − {0}
%   Then, the largest Lyapunov sublevel set Ωc = {x ∈ Rⁿ | V(x) ≤ c}
%   contained in D is an estimate of the region of attraction Rₐ. In other
%   words, once we compute an approximation V(x), we need to check where
%   the approximation is valid, i.e. find D, and then find the largest
%   sublevel set of V(x) in that D. That is your approximation of Rₐ.
%
%   A common practice seems to be to set q(x) = ½ xᵀ Q x even for the
%   nonlinear system [4].
%
%   Summary: Compute an approximation V(x), then check the region where it
%   is valid, then find the largest Lyapunov sublevel set in that region.
%
%   The disconnect between Rₐ and D is what Zubov fixes. In the end, the
%   theorem is very similar to the nonlinear Lyapunov-type equation (2).
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %%%%%%%%%%%%%%%%%%%%% Approach 2: Zubov's Theorem %%%%%%%%%%%%%%%%%%%%%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [Thm 8.10, 5] Let x = 0 be an equilibrium point of (1) and G ⊂ Rⁿ be a
%   domain containing x = 0. Let h:Rⁿ → R be a continuous, positive
%   definite function. Let V:G → R be a continuously differentiable
%   function such that
%        (a) V(0) = 0  and  V(x) > 0  in  G − {0}
%        (b) 0 < V(x) < 1  in  G − {0}
%        (c) As x approaches the boundary of G, lim V(x) = 1
%        (d) 𝜕ᵀV(x)/𝜕x f(x) + h(x)[1 - V(x)] = 0  in  G
%   Then G is the region of attraction Rₐ.
%
%   Notice how in Zubov's theorem, Rₐ = G, rather than dealing with
%   sublevel sets. In other words, instead of having to find a sublevel set
%   that lives in some other domain, we just have to evaluate the 1
%   sublevel set of V(x), and it is not only in G, it is G. So one thing to
%   do is just treat (d) with the Kronecker product approach. This is a
%   significant undertaking, though I am confident I could do it, so I will
%   leave that as a future work for the moment.
%
%   Summary: Not sure how to handle the case where the PDE is not exactly satisfied
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %%%%%%%%%%%%%%% Approach 3: Transformed Zubov's Theorem %%%%%%%%%%%%%%%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   I noticed that in [4], Krener refers to solving the nonlinear Lyapunov
%   equation (2) with q(x) = ½ xᵀ Q x as Zubov's method. Notably, this does
%   not involve the [1 - V(x)] component, and I don't see how this would
%   enforce 0 < V(x) < 1. Digging a little deeper, I found [3] where Tomoki
%   characterizes that they solve 𝜕ᵀV(x)/𝜕x f(x) + h(x)[1 - V(x)] = 0 using
%   polynomial expansions. Therein, the authors mention the substitution
%
%           W(x) = -ln(1 - V(x))    or    V(x) = 1 - exp(-W(x))
%
%   which transforms
%
%           𝜕ᵀV(x)/𝜕x f(x) + h(x)[1 - V(x)] = 0
%           into      𝜕ᵀW(x)/𝜕x f(x) + h(x) = 0
%
%   Notice that the equation for W(x) is the regular nonlinear Lyapunov
%   equation (2). So approach 3 would entail basically solving the same
%   approximation as in approach 1, but instead of playing around with
%   sublevel sets, apply the transformation V(x) = 1 - exp(-W(x)) and just
%   look at the 1 sublevel set. Basically while 0 < V(x) < 1, 0 < W(x) < ထ.
%   Notice how the transformation V(x) = 1 - exp(-W(x)) would map ထ to 1,
%   since  1 - exp(-ထ) = 1. The function exp(-z) should always be a very
%   small nonzero number as z gets very large, and should only equal 0 for
%   z = ထ. However, Matlab numerically returns 1 for 1-exp(-38) due to
%   numerical over/underflow and roundoff type errors. So instead of
%   looking at the 1 sublevel set, we can look at something like the 0.99
%   sublevel set.
%
%   Summary: Compute an approximation W(x), transform it to V(x), check the
%   0.99 sublevel set. As long as this region is in the region where W(x) is
%   a valid Lyapunov function (D), it should be ok.
%
%
%   Authors: Nick Corbin, UCSD
%
%   License: MIT
%
%   Reference: [1] V. I. Zubov, Methods of A.M. Lyapunov and their
%               application. P. Noordhoff, 1964
%              [2] A. M. Lyapunov, “The general problem of the stability of
%               motion," Kharkov Mathematical Society, 1892, doi:
%               10.1080/00207179208934253
%              [3] S. Margolis and W. Vogt, "Control engineering
%               applications of V. I. Zubov's construction procedure for
%               Lyapunov functions," IEEE Transactions on Automatic
%               Control, vol. 8, no. 2, pp. 104–113, Apr. 1963, doi:
%               10.1109/tac.1963.1105553
%              [4] A. Krener, C. Aguilar, and T. Hunt, “Series solutions of
%               HJB equations,” Mathematical System Theory—Festschrift in
%               Honor of Uwe Helmke on the Occasion of his Sixtieth
%               Birthday, pp. 247–260, 2013
%              [5] H. K. Khalil, Nonlinear systems, Third edition, Pearson
%               Education, 2013
%
%  Part of the PPR repository.
%%

if (nargin < 4)
    verbose = false;
    if (nargin < 3)
        degree = length(f);
    end
end

n = length(f{1});

%% Method 1: Solve nonlinear Lyapunov equation for V(x), then find largest sublevel set in D
% Solve 𝜕ᵀV(x)/𝜕x f(x) + q(x) = 0 using ppr()
options.verbose = verbose;
[v] = ppr(f, zeros(n,n), q, 1/inf, degree, options);

% TODO: Find largest sublevel set in D

%% Method 2: Solve proper Zubov equation for V(x), then find 1 sublevel set
% TODO

%% Method 3: Solve nonlinear Lyapunov equation for W(x), transform, then find 1 sublevel set
% Solve 𝜕ᵀV(x)/𝜕x f(x) + q(x) = 0 using ppr()
options.verbose = verbose;
[v] = ppr(f, zeros(n,n), q, 1/inf, degree, options);

% TODO: Find largest sublevel set in D

end
