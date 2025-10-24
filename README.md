# PPR Examples
Polynomial and analytic examples for testing Polynomial-Polynomial Regulator (PPR) design. 

## Syntax
`[f, g, h] = getSystemXX()`

## Description
The PPR problem computes polynomial approximations to the value function and optimal feedback law for nonlinear systems featuring analytic or polynomial nonlinearities. 
The syntax for PPR is designed to be similar to Matlab's `lqr()`, but instead of passing system matrices $A$ and $B$, we pass $f(x)$ and $g(x)$.
The PPR code is written to handle these functions in polynomial form, so instead of passing matrices `A` and `B`, we pass cell arrays `f` and `g` containing the polynomial coefficients for $f(x)$ and $g(x)$.

In this directory, we have many different examples of common benchmark problems for nonlinear control design. 
Each example features a function `[f, g, h] = getSystemXX()` to return the polynomial coefficients for the drift and input functions. 
There is also generally a third output for an output equation, but the PPR problem is a state feedback problem so this is not generally used. 
Below, we briefly summarize the models contained in each example.

## Examples
- System 1: An academic 1D example for which the nonlinear balanced truncation energy functions are known analytically. 
- System 2: An academic 2D quadratic bilinear model, also used for computing nonlinear balancing energy functions. 
- System 3: Finite element model of Burger's equation in 1D.
- System 4: Finite element model of Kuramoto-Sivashinky equation in 1D.
- System 5: Unicycle quadratic-bilinear approximation.
- System 6: Finite element nonlinear Euler-Bernoulli beam.
- System 7: Aircraft stall stabilization example.
- System 8: Nonlinear reaction-diffusion heat equation.
- System 9: 1D Allen-Cahn equation, pseudospectral Chebychev discretization. 
- System 10: An academic 1D example with transcritical bifurcation.
- System 11: Inverted pendulum model (not nondimensionalized; see System 26).
- System 12: An academic 2D model cooked up for balancing.  
- System 13: An academic 2D model cooked up for balancing.  
- System 14: 2D gradient system describing a double pendulum.
- System 15: 4D double pendulum model.
- System 16: 6D triple pendulum model.
- System 17: Coupled chain of mass-spring-dampers
- System 18: A linear 3D academic model.
- System 19: Deprecated, since deleted.
- System 20: Vlasov equation.
- System 21: 1D Allen-Cahn equation, pseudospectral Chebychev discretization. 
- System 22: 4D inverted pendulum on a cart. 
- System 23: An academic 3D model for model reduction.
- System 24: A 2D academic model, simplified version of System 23.
- System 25: Van der Pol oscillator.
- System 26: Nondimensionalized 2D inverted pendulum model.
- System 27: 1D unstable heat equation finite element model.
- System 28: 2D LINEAR heat equation model.
- System 29: 2D Allen-Cahn equation polynomial finite element model. Optimized for sparsity. 
- System 30: 1D nonlinear finite element heat equation model. 







