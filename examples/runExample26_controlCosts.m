function runExample26_controlCosts()
%runExample26_controlCosts Runs 2D inverted pendulum example (normalized).
%   Computes the control costs for the table.
%
%   Usage:  runExample26_controlCosts()
%
%   Description: The basic equation for the inverted pendulum is
%
%               ẍ - sin(x) = u(t)
%
%   This can be put in first-order form with x₁ = x, x₂ = ẋ as
%
%             ẋ₁ = x₂
%      (1)    ẋ₂ = sin(x₁) + u(t)
%              y = x₁
%
%   Now, we can compute a controller u(x) = K(x) using PPR.
%
%   References: [1]
%
%   Part of the PPR repository.
%%
fprintf('Running Example 26\n')

% Get dynamics
degree = 12;
[f, g, ~] = getSystem26(degree-1,false);
F = @(x) [x(2); sin(x(1))]; G = @(x) g{1};

% Get value function/controller
q = 0; R = 1;
[v, K] = ppr(f, g, q, R, degree);
% [v, K] = ppr(f, g, q, R, degree, options);

%% Simulate and compute control costs
fprintf('\n\n# Table 1 Data (Pendulum Control Costs)\n');
fprintf('# Control costs for different controllers\n');
fprintf("    d  &      V(x)    &  Integrated Cost     \n")

x0 = [-pi;2]; tspan = [0 10];
for d=2:2:degree
    [t, X] = ode45(@(t, x) F(x) + G(x) * kronPolyEval(K, x, degree=d-1), tspan, x0);
    
    % Compute performance index (cost)
    runningCost = zeros(size(t));
    for i=1:length(t)
        x = X(i,:).'; Ux = kronPolyEval(K, x, degree=d-1);
        runningCost(i) = 1/2*(Ux.'*R*Ux);
    end
    integratedCost = trapz(t, runningCost);
    
    fprintf("   %2i  &   %8.5f   &      %8.5f    \n", d, 0.5*kronPolyEval(v, x0, degree=d), integratedCost)
end
fprintf("\n")


end



