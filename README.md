# ppr
Polynomial-Polynomial Regulator (PPR) design 

## Syntax
`[v,K] = ppr(f,g,q,r,d)`

`[v,K,options] = ppr(f,g,q,r,degree,options)`

## Description
`[v,K] = ppr(f, g, q, r, d)` calculates the optimal gain coefficients `K` and the first `d` value function coefficients `v` of the associated Hamilton-Jacobi-Bellman equation using the continuous-time state-space model defined by `f` and `g`. `q` and `r` are the weight coefficients for states and inputs, respectively. 
All of these quantities are cell arrays containing matrix coefficients defining polynomials (feedback law, value function, dynamics, cost function, etc.).

<!-- example -->

## Examples
<!-- collapse all -->
<!-- ### LQR Control for Inverted Pendulum Model -->

<details open>
<summary>PPR Control for Inverted Pendulum on a Cart Model</summary>
<br>

`getSystem22()` returns the Taylor approximation of the 4D state-space model of an inverted pendulum on a cart based on [5-7].
The outputs are the cart displacement $x$ and the pendulum angle $\theta$. 
The control input $u$ is the horizontal force on the cart.
The nonlinear equations of motion are

```math
\begin{bmatrix}\dot{x}\\\\ddot{x}\\\\dot{\theta}\\\\ddot{\theta}\end{bmatrix}=\begin{bmatrix}x_2\\\\frac{-0.27x_2-0.18x_4^2\sin(x_3)+4\cos(x_3)\sin(x_3)}{1.8-0.44\cos(x_3)^2}\\\x_4\\\\frac{40\sin(x_3)-0.67x_2\cos(x_3)-0.44x_4^2\cos(x_3)\sin(x_3)}{1.8-0.44\cos(x_3)^2}\end{bmatrix}+\begin{bmatrix}0\\\\frac{2.7}{1.8-0.44\cos(x_3)^2}\\\0\\\\frac{6.7\cos(x_3)}{1.8-0.44\cos(x_3)^2}\end{bmatrix}u
```

```math
y =
    \begin{bmatrix}
        1 & 0 & 0 & 0 \\
        0 & 0 & 1 & 0
    \end{bmatrix}
    \begin{bmatrix}
        x \\
        \dot{x} \\
        \theta \\
        \dot{\theta}
    \end{bmatrix}
    +
    \begin{bmatrix}
        0 \\
        0
    \end{bmatrix} u
```

The linearized dynamics are then
```math
\begin{bmatrix}
        \dot{x} \\
        \ddot{x} \\
        \dot{\theta} \\
        \ddot{\theta}
    \end{bmatrix}
    =
    \begin{bmatrix}
        0 & 1 & 0 & 0 \\
        0 & -0.2 & 3 & 0 \\
        0 & 0 & 0 & 1 \\
        0 & -0.5 & 30 & 0
    \end{bmatrix}
    \begin{bmatrix}
        x \\
        \dot{x} \\
        \theta \\
        \dot{\theta}
    \end{bmatrix}
    +
    \begin{bmatrix}
        0 \\
        2 \\
        0 \\
        5
    \end{bmatrix} u
```
which is almost identical to the model used in the LQR function documentation, except for the (2,2) entry which is -0.1 there. 
`getSystem22()` uses the helper function `approxPolynomialDynamics()` to compute the polynomial coefficient arrays `f`, `g`, and `h` for the dynamics using the symbolic nonlinear dynamics.

Run the example with `runExample22()` after running `setKroneckerToolsPath`, `addpath('examples')`, and `addpath('utils')`. Here is a breakdown of key steps in the example script:

First, load a degree 7 approximation to the state-space model to the workspace. `f` and `g` are cell arrays containing the Taylor approximations to the nonlinear dynamics, and `xdot` is a symbolic function capturing the full nonlinear dynamics for use in simulations.
```
[f, g, ~, xdot] = getSystem22(7);
```

Define the initial condition of interest.
Neither PPR nor LQR can do full swing-up control currently, but there are conditions where PPR succeeds while LQR fails.
A simple case is when the cart is already in the up position, but shifted to the right 10 units and we wish to move it to the left to the origin without tipping the pendulum over.
```
x0 = [12;0;0;0];
```
Another interesting case is `x0 = [0;0;-pi/3*1.04;0];`, which physically corresponds to the pendulum starting from a tilted position. 
This case is just beyond the region where LQR succeeds; however, PPR still succeeds at this condition.

Next, we can define and compute the control laws.
Let's use a simple quadratic cost function defined by $Q$ and $R$ matrices: 
```
Q = [1,0,0,0;...
    0,0,0,0;...
    0,0,1,0;...
    0,0,0,0];
R = 1;
```
Find the gain matrices $K_i$ using `ppr`. 

```
[~, K] = ppr(f, g, Q, R,8)
```
    K = 1×7 cell array

        {1×4 double}    {1×16 double}    {1×64 double}    {1×256 double}    {1×1024 double}    {1×4096 double}    {1×16384 double}


Unlike `lqr`, which returns a single gain matrix $K$, `ppr` return several matrices in a cell array. 
The first entry in the array $K_1$ is precisely the 1x4 LQR gain matrix, hence an LQR controller can be computed by calling `ppr` with `degree` set to 2, which computes a quadratic value function approximation and a linear feedback. 
The remaining $K_i$ are the higher order gain matrices corresponding to quadratic, cubic, etc. feedback terms. 

A feedback law can be defined using the gain coefficient cell array `K` using an anonymous function based on the utility function `kronPolyEval'. 
Here we use the function to define both the LQR and PPR controllers, along with an "open-loop" controller corresponding to the uncontrolled system: 
```
uOpenLoop = @(z) zeros(1,1);
uLQR = @(z) (kronPolyEval(K, z, 1));
uPPR = @(z) (kronPolyEval(K, z));
```

Using this format, we can easily simulate the closed-loop performance of the different controllers. 
Define the system dynamics $\dot{x} = f(x,u)$ using an anonymous function based on the `xdot` function from earlier:
```
FofXU = @(x,u) xdot(x,u);
```

The closed-loop systems can be simulated with your favorite ode solver; for example, 
```
% Time vector
tmax = 10; dt = 0.001; t = 0:dt:tmax-dt;

% Simulate with ode45
[t, Xunc] = ode45(@(t, x) FofXU(x,uOpenLoop(x)),t, x0);
[t, XLQR] = ode45(@(t, x) FofXU(x,uLQR(x)),t, x0);
[t, XPPR] = ode45(@(t, x) FofXU(x,uPPR(x)),t, x0);

% Plot the results 
figure; 
for i=1:4; subplot(2,2,i); hold on; plot(t,Xunc(:,i)); end
for i=1:4; subplot(2,2,i); hold on; plot(t,XLQR(:,i)); end
for i=1:4; subplot(2,2,i); hold on; plot(t,XPPR(:,i)); end

subplot(2,2,1); xlabel('time'); ylabel('z'); title('cart position'); ylim([-2 15])
subplot(2,2,2); xlabel('time'); ylabel('zdot'); title('cart velocity'); ylim([-10 5])
subplot(2,2,3); xlabel('time'); ylabel('theta'); title('pendulum angle'); ylim([-2*pi 2*pi])
subplot(2,2,4); xlabel('time'); ylabel('thetadot'); title('pendulum velocity'); ylim([-3 4])
legend('Uncontrolled system','LQR','PPR','Location','southoutside')
```

`runExample22()` does essentially this, but also compares with nonlinear MPC based on Matlab documentation's example.
To make the comparison easier, an explicit Euler time integrator is used, which slightly changes the results, but they are qualitatively similar. 
Below is the plot of the results comparing the controller performances: 

<img src="https://cnick1.github.io/images/invertedPendulumExample.png" width="600">
<img src="https://cnick1.github.io/images/invertedPendulumExampleControls.png" width="400">

This is not necessarily the optimal MPC performance that can be achieved with modern solvers; this is just based on modifying the simple example in the Matlab documentation. 
It is worth noting that the controller costs here are 180.642 for PPR and 204.601 for MPC, so PPR performs better. 
However, PPR is significantly faster online, because the control law is given as a simple polynomial, whereas MPC must solve an optimization problem repeatedly. 
Better performance can of course be achieved by modifying the settings used for MPC, but then the computations become more expensive and take even longer. 
It is also worth noting that MPC can do swing-up control, whereas PPR and LQR currently fail.

</details>

<details open>
<summary>PPR Control for Aircraft Stabilization Model</summary>
<br>
`getSystem7()` returns the cubic 3D state-space model of an F-8 aircraft cruising at 30,000 ft at Mach = 0.85 developed originally in [5].
The state represents the angle of attack $x_1$, the angle f the plane relative to the trim pitch $x_2$, and the rotation rate of the plane $x_3$.
The control input $u$ is the angle of the tail elevator.
The nonlinear equations of motion are

$$
    \dot{x}_1  = x_3  - x_1^2 x_3 - 0.088 x_1 x_3 - 0.877 x_1 + 0.47 x_1^2 - 0.019 x_2^2 + 3.846 x_1^3 - 0.215 u + 0.28 u x_1^2\\
    \dot{x}_2  = x_3  \\                                           
    \dot{x}_3  = -0.396 x_3 - 4.208 x_1 - 0.47 x_1^2 - 3.564 x_1^3 - 20.967 u + 6.265 u x_1^2.
$$ 

The objective is to stabilize the aircraft subject to a perturbation in the initial angle of attack, i.e. bring the initial condition $x_0 = [\alpha_0 \quad 0 \quad 0]^\top$ asymptotically to the origin.

Run the example with `runExample7()` after running `setKroneckerToolsPath` and `addpath('examples')`. Here is a breakdown of key steps in the example script:

First, load the state-space model to the workspace. `f` and `g` are cell arrays containing the coefficients to the polynomial dynamics.
```
[f, g, ~] = getSystem7();
```

Define the initial condition of interest. 
Interesting choices are $\alpha_0 = 25, 27, 30, 35$.

```
alpha0 = 27; 
x0 = [pi / 180 * alpha0; 0; 0];
```

Next, we can define and compute the control laws.
Let's use a simple quadratic cost function defined by $Q$ and $R$: 
```
Q = 0.25; R = 1; 
```
Note that if a scalar is passed as either $Q$ or $R$, it will be converted to a diagonal matrix automatically, hence $Q$ here is really a 3 by 3 matrix.

Now compute a degree 7 approximation to the optimal control using `ppr`; degree 1, 3, and 5 controllers are given by subsets of the gain matrices $K_i$.
```
[~, K] = ppr(f, g, Q, R,8)
```
    K = 1×7 cell array

        {1×3 double}    {1×9 double}    {1×27 double}    {1×81 double}    {1×243 double}    {1×729 double}    {1×2187 double}


Unlike `lqr`, which returns a single gain matrix $K$, `ppr` return several matrices in a cell array. 
The first entry in the array $K_1$ is precisely the 1x3 LQR gain matrix, hence an LQR controller can be computed by calling `ppr` with `degree` set to 2, which computes a quadratic value function approximation and a linear feedback. 
The remaining $K_i$ are the higher order gain matrices corresponding to quadratic, cubic, etc. feedback terms. 

A feedback law can be defined using the gain coefficient cell array `K` using an anonymous function based on the utility function `kronPolyEval'. 
Here we use the function to define both the LQR and PPR controllers, along with an "open-loop" controller corresponding to the uncontrolled system: 
```
uLQR = @(z) (kronPolyEval(K, z, 1));
uPPR3 = @(z) (kronPolyEval(K, z, 3));
uPPR5 = @(z) (kronPolyEval(K, z, 5));
uPPR7 = @(z) (kronPolyEval(K, z, 7));
```

Using this format, we can easily simulate the closed-loop performance of the different controllers. 
Define the system dynamics $\dot{x} = f(x) + g(x)u $ using anonymous functions:
```
F = @(x) kronPolyEval(f, x);
G = @(x) (g{1} + kronPolyEval(g(2:end), x));
```

The closed-loop systems can be simulated with your favorite ode solver; for example, 
```
% Time vector
tspan = [0, 12]; 

% Simulate with ode45
[tLQR, XLQR]   = ode45(@(t, x) F(x) + G(x) * uLQR(x),      tspan, x0);
[tPPR3, XPPR3] = ode45(@(t, x) F(x) + G(x) * uPPR3(x),     tspan, x0);
[tPPR5, XPPR5] = ode45(@(t, x) F(x) + G(x) * uPPR5(x),     tspan, x0);
[tPPR7, XPPR7] = ode45(@(t, x) F(x) + G(x) * uPPR7(x),     tspan, x0);

% Plot the results 
figure; hold on;
plot(tLQR,  XLQR(:,1)  / (pi/180))
plot(tPPR3, XPPR3(:,1) / (pi/180))
plot(tPPR5, XPPR5(:,1) / (pi/180))
plot(tPPR7, XPPR7(:,1) / (pi/180))

xlabel('time'); ylabel('\alpha, degrees'); title('Aircraft angle of attack'); ylim([0 40])
legend('LQR','Degee 3 PPR','Degree 5 PPR','Degree 7 PPR')
```

Below is the plot of the results comparing the controller performances for different initial angles of attack: 

<img src="https://cnick1.github.io/images/aircraftStallStabilizationExample.png" width="1000">

</details>

## Input Arguments
<!-- collapse all -->
<details open>
<summary>f — Drift vector field coefficients</summary>
<br>
cell array

Drift vector field, specified as a cell array containing the polynomial coefficients of the drift in Kronecker polynomial form. 
The first coefficient is the n-by-n linear state matrix, `A = f{1}`, where n is the number of states. 
The remaining coefficients `f{i}` must be n-by-n^i matrices.
</details>

<details open>
<summary>g — Input vector field coefficients</summary>
<br>
cell array

Input vector field(s), specified as a cell array containing the polynomial coefficients of the input-to-state map in Kronecker polynomial form. 
The first coefficient is the n-by-m linear input-to-state matrix, `B = g{1}`, where m is the number of inputs.
The remaining coefficients `g{i}` must be n-by-mn^(i-1) matrices.
</details>

<details open>
<summary>q — State-cost coefficients</summary>
<br>
cell array

State-cost coefficients, specified as a cell array containing the polynomial coefficients of the state penalty term in the cost function in Kronecker polynomial form.
The first entry is ignored and should be empty, `q{1} = []`.
The second entry is the quadratic state penalty, typically defined as an n-by-n matrix `q{2} = Q`, where n is the number of states.
Alternatively a scalar can be passed, in which case it is automatically multiplied by an n-by-n identity matrix.
The remaining coefficients `q{i}` must be n^i-by-1 vectors. 
Alternatively scalars can be passed, in which case they are automatically multiplied by n-by-n-by-...-by-n identity tensors.
The quadratic cost term `q{2}` can also be vectorized for consistency.

Optionally, just the quadratic cost matrix `Q` can be passed.

`q{2} = Q` should be a positive semi-definite matrix.
</details>

<details open>
<summary>r — Input-cost coefficients</summary>
<br>
cell array 

Input-cost coefficients, specified as a cell array containing the polynomial coefficients of the input penalty term in the cost function in Kronecker polynomial form.
The first entry is the constant-in-state input penalty, typically defined as an m-by-m matrix `r{1} = R`, where m is the number of inputs.
Alternatively a scalar can be passed, in which case it is automatically multiplied by an m-by-m identity matrix.
The remaining coefficients `r{i}` must be m^2-by-n^(i-1) matrices. 
Alternatively scalars can be passed, in which case they are automatically multiplied by m-by-m-by-...-by-m identity tensors.
The first cost term `r{1}` can also be vectorized for consistency.

Optionally, just the quadratic cost matrix `R` can be passed.

`r{1} = R` should be a positive definite matrix.
</details>

<details open>
<summary>d — Value function polynomial degree</summary>
<br>
length(f)+1 (default)| (typically even) integer

Value function desired polynomial degree, typically specified as an even integer.
A degree `d` value function corresponds to a degree `d-1` approximation to the optimal feedback law. 
A degree `d` value function approximation only includes terms from the cost function of degree `d` or less and the dynamics of degree `d-1` or less, so you don't have to provide infinitely many terms in the `f`, `g` cell arrays.
The default choice of `d` is `length(f)+1`.
</details>


<details open>
<summary>options — Solver options</summary>
<br>
structure

Solver options with the following fields:
- `verbose` - boolean determining if runtime information is printed.
- `r` - dimension of the reduced-order model used for accelerating the computation of higher-order coefficients. (Defaults to $n$, the full state dimension.) This option can be used to accelerate PPR, although the computed solutions will no longer be true Taylor coefficients of the value function and optimal control law. The reduction uses the $r$ leading eigenvectors of the quadratic value function coefficient $V_2$ coming from the LQR problem on the linearized dynamics as a reduction basis. The intuition behind this is that small eigenvalues in $V_2$ correspond to directions in which the gradient of the value function is locally small, meaning that those directions contribute little to the controller's performance. Hence those directions can be omitted, and the higher-order feedback can be computed only in the directions that require large control efforts for the linearized system.
- `method` - model reduction method for accelerating the computation of higher-order coefficients. The default is to use the dominant eigenvectors of $V_2$.
- `fr`,`gr` - reduced-order dynamics of dimension $r$ if they have already been computed.

</details>




## Output Arguments
<!-- collapse all -->
<details open>
<summary>v — Coefficients of solution approximation to the associated Hamilton-Jacobi-Bellman equation</summary>
<br>
cell array

Taylor coefficients of the value function, which is solution to the associated Hamilton-Jacobi-Bellman equation, specified in Kronecker polynomial form.
`v{1}` is empty; `v{2}=V2` is the solution to the algebraic Riccati equation associated with the LQR problem on the linearized dynamics. 
`v{i}' is an n^i-by-1 vector.
</details>

<details open>
<summary>K — Optimal gain coefficients</summary>
<br>
cell array 

Taylor coefficients of the optimal feedback law, specified in Kronecker polynomial form. 
`K{1}` is the m-by-n LQR gain matrix, where m is the number of inputs and n is the number of states.
The higher-order `K{i}`, which are m-by-n^i matrices, are the higher-order coefficients in the Taylor expansion of the optimal control law.
</details>

<details open>
<summary>options — Solver options</summary>
<br>
structure

Solver options and information with the following fields:
- `verbose` - boolean determining if runtime information is printed.
- `r` - dimension of the reduced-order model used for accelerating the computation of higher-order coefficients. 
- `fr`,`gr` - reduced-order dynamics of dimension $r$ if they have been computed so that they can be reused if necessary.
- `Tib`,`TibInv` - transformation to and from reduced-order state

</details>



## Limitations
The input data must satisfy the following conditions:
- The pair (`f{1}`,`g{1}`) must be stabilizable.
- $R$ must be positive definite.
- $Q$ must be positive semidefinite ($Q$≥0).
These assumptions essentially correspond to requiring the LQR problem on the linearized dynamics to have a solution. 

The solution approximations computed herein are of a Taylor-series nature; hence, they are guaranteed to exist if the above assumptions are satisfied, and they are guaranteed to converge locally. 
However, the region of convergence is not guaranteed to be global, so the solution may only be valid in a neighborhood of the origin. 
In practice, the more strongly nonlinear the problem is, the more localized you can expect the solution to be.

Since the solutions computed herein are of a polynomial nature, the value function is expected to blow up to negative or positive infinity. 
The LQR value function which comes from the first term in the Taylor expansion and is a quadratic function guaranteed to be positive, satisfying the requirements of a Lyapunov function. 
The higher-order polynomial approximations to the value function however are not guaranteed to satisfy the positivity requirement of a Lyapunov function outside of a neighborhood of the origin. 

## Algorithms
For continuous-time systems, `ppr` computes the multi-input state-feedback control $u=K_1x+K_2x^{(2)}+K_3x^{(3)}$---where
 $x^{(2)}=(x\otimes x),x^{(3)}(x\otimes x\otimes x),etc.$ and $\otimes$ is the Kronecker product---that minimizes the polynomial cost function
 
$$\frac{1}{2} \int_{0}^{\infty} \left(x^\top Q x  +  u^\top R u + \sum_{p=3}^{\infty} q_p^\top x^{(p)} + \frac{1}{2} \sum_{p=1}^{\infty}
  {x^{(p)}}^\top r_p^\top u^{(2)}\right)dt$$
  
subject to the system dynamics 

$$\dot{x}  = \sum_{p=1}^\infty F_p x^{(p)} + \left(\sum_{p=0}^\infty G_p \left(x^{(p)} \otimes I_m\right) \right) u.$$

This can be thought of as a Taylor approximation to the analytic optimal control problem

$$J(x,u)=\int_0^\infty \left(x^\top Q(x) x + u^\top R(x) u \right)dt$$

subject to the control-affine system dynamics 

$$\dot{x} = f(x)+ g(x)u.$$

In addition to the state-feedback gains $K_i$, ppr returns the Taylor coefficients $v_i$ of value function associated with the Hamilton-Jacobi-Bellman equation

$$0 = \frac{\partial V^\top(x)}{\partial x}  f(x) - \frac{1}{2} \frac{\partial V^\top(x)}{\partial x}g(x) R^{-1}(x)g^\top(x) \frac{\partial V(x)}{\partial x} + \frac{1}{2} x^\top Q(x) x.$$

The gain coefficients $K_i$ are derived from $v_i$ using

$$u(x)  =  -R^{-1}(x) g^\top(x) \frac{\partial V(x(t))}{\partial x}.$$


## Dependencies
The PPR repository is dependent on the `KroneckerTools` repository and uses some utility functions from the `tensor_toolbox` repository.
Clone these repositories into the same parent folder as PPR, or modify the path in `setKroneckerToolsPath.m`.
```
  git clone https://www.github.com/cnick1/PPR.git
  git clone https://www.github.com/cnick1/KroneckerTools.git
  git clone https://gitlab.com/tensors/tensor_toolbox.git
```

`ppr` uses a structured linear solver when computing the higher-order coefficients. 
Either download the `tensor_recursive` package from https://www.epfl.ch/labs/anchp/index-html/software/misc/ to use the more efficient solver, or change the solver from `solver = 'chen-kressner'` to `solver = 'bartels-stewart'` in `KroneckerSumSolver.m` in the `KroneckerTools` repository. 

## References 
[1] N. A. Corbin and B. Kramer, “Computing solutions to the polynomial-polynomial regulator problem,” in 2024 63rd IEEE Conference on Decision and Control, Dec. 2024. doi: 10.48550/arXiv.2410.22291

[2] J. Borggaard and L. Zietsman, “The quadratic-quadratic regulator problem: approximating feedback controls for quadratic-in-state nonlinear systems,” in 2020 American Control Conference (ACC), Jul. 2020, pp. 818–823. doi: 10.23919/ACC45564.2020.9147286

[3] J. Borggaard and L. Zietsman, “On approximating polynomial-quadratic regulator problems,” IFAC-PapersOnLine, vol. 54, no. 9, pp. 329–334, 2021, doi: 10.1016/j.ifacol.2021.06.090

[4] M. Chen and D. Kressner, “Recursive blocked algorithms for linear systems with Kronecker product structure,” Numerical Algorithms, vol. 84, no. 3, pp. 1199–1216, Sep. 2019, doi: 10.1007/s11075-019-00797-5.

[5] https://www.mathworks.com/help/control/ref/lti.lqr.html

[6] https://www.mathworks.com/help/symbolic/derive-and-simulate-cart-pole-system.html

[7] https://www.mathworks.com/help/mpc/ug/swing-up-control-of-a-pendulum-using-nonlinear-model-predictive-control.html

[8] W. L. Garrard and J. M. Jordan, “Design of nonlinear automatic flight control systems,” Automatica, vol. 13, no. 5, pp. 497–505, Sep. 1977, doi: 10.1016/0005-1098(77)90070-x

## Version History
- v0.9.1 - Updated release associated with v1.0.0 of cnick1/NLbalancing for IEEE TAC.
- v0.9.0 - Initial version of examples submitted to CDC2024. 

## Credits
Developed and maintained by Nicholas Corbin. 
- Inspired heavily by Jeff Borggaard's QQR and PQR work and associated repositories, but largely rewritten. 
`cnick1/KroneckerTools` is a fork of Jeff Borggaard's `KroneckerTools` repository.
- Credit to Minhong Chen and Daniel Kressner for a wonderful Lyapunov structure linear solver.
- Credit to Brett W. Bader, Tamara G. Kolda, Daniel M. Dunlavy, et al. for the spectacular tensor_toolbox. 

## See Also
`kronPolyEval` | `kronPolyDerivEval` 
