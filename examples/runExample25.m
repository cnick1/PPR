function runExample25(lienard)
%runExample25 Runs 2D Van der Pol example to compare with [2] and [3].
%
%   Usage:  runExample25()
%
%   Inputs:    lienard - whether to use model (1) (false) or (2) (true)
%        
%   Description: The Van der Pol oscillator is a nonlinear oscillator with
%   nonlinear damping. The basic equation in second-order form is 
%
%               ẍ + ϵ(1 - x²)ẋ + x = u(t)
%
%   This can be put in first-order form with x₂ = ẋ as   
%
%             ẋ₁ = x₂
%      (1)    ẋ₂ = -x₁ - ϵ(1 - x₁²)x₂ + u(t)
%              y = x₁
% 
%   Alternatively, [3] defines x₂ = ẋ + ϵ(x - x³/3) to obtain
% 
%             ẋ₁ = x₂ - ϵ(x₁ - x₁³/3)
%      (2)    ẋ₂ = -x₁
%              y = x₁
% 
%   For lienard = true, the results to compare with [3] are produced. For
%   lienard = false, the values from [2] are used; however, note that in
%   that case the plots are mirrored compared to [2], which technically
%   considers the reverse-time unstable Van der Pol oscillator rather than
%   the stable case. 
%
%   References: [1] https://en.wikipedia.org/wiki/Van_der_Pol_oscillator
%               [2] H. K. Khalil, Nonlinear systems, Third edition, Pearson
%                Education, 2013
%               [3] S. Margolis and W. Vogt, "Control engineering
%                applications of V. I. Zubov's construction procedure for
%                Lyapunov functions," IEEE Transactions on Automatic
%                Control, vol. 8, no. 2, pp. 104–113, Apr. 1963, doi:
%                10.1109/tac.1963.1105553
%
%   Part of the PPR repository.
%%
if nargin < 1 
    lienard = false;
end

fprintf('Running Example 25\n')

if lienard
    clin = 0.875;
    % cnonlin = 0.985; degree = 6;
    cnonlin = 0.995; degree = 10;
    Q = [1 0; 0 0];
    eps = 0.7;
else
    % clin = 0.9;
    % cnonlin = 0.9975;
    % Q = [1 0; 0 1];
    % eps = 1;

    clin = 0.465;
    cnonlin = 0.975;
    Q = [1 0; 0 0];
    eps = 0.7;
    degree = 10;
end

[f, g, h] = getSystem25(lienard, eps);

% Using PPR 
v = zubov(f, Q, degree);

nX = 101; nY = nX; xLim = 3; yLim = 3;
xPlot = linspace(-xLim, xLim, nX); yPlot = linspace(-yLim, yLim, nY); [X, Y] = meshgrid(xPlot, yPlot);

Vlin = zeros(nY, nX); 
Vdotlin = zeros(nY, nX); 
Vnonlin = zeros(nY, nX); 
Vdotnonlin = zeros(nY, nX); 

for i = 1:nY
    for j = 1:nX
        x = [X(i, j); Y(i, j)];
        
        Vlin(i,j) = kronPolyEval(v,x,degree=2);
        Vdotlin(i,j) = kronPolyDerivEval(v,x,2) * kronPolyEval(f,x);


        Vnonlin(i,j) = kronPolyEval(v,x);
        Vdotnonlin(i,j) = kronPolyDerivEval(v,x) * kronPolyEval(f,x);

    end
end

cmap = [1, 0.6, 0.6;  % Light red (V < 0)
        0.6, 1, 0.6]; % Light green (V > 0)

%% Plot region of attraction using quadratic Lyap estimate 
figure; hold on; xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex'); title('R_A estimate using quadratic Lyapunov function')
% Define colors: light green for V > 0, red for V < 0
colormap(cmap)

% Plot region where V is positive definite 
[c,h] = contourf(X, Y, Vlin,[-inf 0]); h.FaceAlpha = 0.25;

% Plot region where Vdot is negative definite 
[c,h] = contourf(X, Y, -Vdotlin,[-inf 0]); h.FaceAlpha = 0.25;

% Plot largest Lyapunov level set that fits in the region
ax1 = gca; xlim([-3 3]); ylim([-3 3]);
ax2 = axes; xlim([-3 3]); ylim([-3 3]);

contour(X, Y, 1-exp(-Vlin), clin*[1 1],'k'); ax2.Visible = 'off';

plotLimitCycle(@(x) kronPolyEval(f,x), [-.1,.1], [-.1,.1], [0,100], 2)

%% Plot region of attraction using nonlinear Lyap estimate 
figure; hold on; xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex'); title('R_A estimate using polynomial Lyapunov function')
% Define colors: light green for V > 0, red for V < 0
colormap(cmap)

% Plot region where V is positive definite 
[c,h] = contourf(X, Y, Vnonlin,[-inf 0]); h.FaceAlpha = 0.25;

% Plot region where Vdot is negative definite 
[c,h] = contourf(X, Y, -Vdotnonlin,[-inf 0]); h.FaceAlpha = 0.25;

% Plot largest Lyapunov level set that fits in the region 
ax1 = gca; xlim([-3 3]); ylim([-3 3]);
ax2 = axes; xlim([-3 3]); ylim([-3 3]);

contour(X, Y, 1-exp(-Vnonlin), cnonlin*[1 1],'k'); ax2.Visible = 'off';    

plotLimitCycle(@(x) kronPolyEval(f,x), [-.1,.1], [-.1,.1], [0,100], 2)

%% Plot Ra approximations
figure; hold on; xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex'); title('R_A estimates')
contour(X, Y, 1-exp(-Vlin), clin*[1 1],'k:'); 
contour(X, Y, 1-exp(-Vnonlin), cnonlin*[1 1],'k');     

plotLimitCycle(@(x) kronPolyEval(f,x), [-.1,.1], [-.1,.1], [0,100], 2)

if lienard
    xlim([-2.4 4]); ylim([-2.2 2.5]);
else
    axis equal; xlim([-3 3]); ylim([-3 3]);
    legend('R_A d=2 approximation','R_A d=10 approximation','Limit cycle')
end


end

%% Helper functions
function plotLimitCycle(f, x_range, y_range, t_span, num_points)
% plot_phase_portrait - Plot limit cycle by integrating backwards in time
%
% Inputs:
%   f        - Function handle representing dx/dt = f(x), where x is a 2D vector
%   x_range  - [xmin, xmax] range for initial x values
%   y_range  - [ymin, ymax] range for initial y values
%   t_span   - Time span for integration (e.g., [0, 10])
%   num_points - Number of initial conditions in each direction

[X0, Y0] = meshgrid(linspace(x_range(1), x_range(2), num_points), ...
    linspace(y_range(1), y_range(2), num_points));

hold on;
for i = 1:numel(X0)
    x0 = [X0(i); Y0(i)];
    [~, Xb] = ode45(@(t, x) -f(x), t_span, x0);
    plot(Xb(end-100:end, 1), Xb(end-100:end, 2), 'r');
end

end
