function runExample26(degree)
%runExample26 Runs 2D inverted pendulum example (normalized). Plots a) the value
%   function, b) the HJB residual, and c) closed-loop phase portraits.
%
%   Usage:  runExample26()
%
%   Inputs:    degree - degree of value function approximation
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
%   and sin(x₁) can be approximated as
%
%             sin(x₁) = x₁ - x₁³/6 + x₁⁵/120 - x₁⁷/5040 + x₁⁹/362880 + ...
%
%   In this script, we compute the value function, the HJB residual, and
%   closed-loop phase portraits to understand the nature of the polynomial
%   approximations in the context of optimal control.
%
%   References: [1]
%
%   Part of the PPR repository.
%%
if nargin < 1
    degree = 6;
end

fprintf('Running Example 26\n')

[f, g, ~] = getSystem26(degree-1);

F = @(x) [x(2); sin(x(1))];
G = @(x) g{1};

%% Open-loop phase portrait
plotPhasePortrait(1, @(x) 0)

%% Value function & HJB Residual plots
%  Compute the polynomial approximations to the value function
q = 0; r = 1; tic;
[v, K, ~] = ppr(f, g, q, r, degree);
fprintf("Completed ppr() in %2.2f seconds. \n", toc)

nX = 301; nY = nX; xLim = pi; yLim = 5;
xPlot = linspace(-xLim, xLim, nX); yPlot = linspace(-yLim, yLim, nY); [X, Y] = meshgrid(xPlot, yPlot);

HJBResidual = zeros(nY, nX); valueFunction = zeros(nY, nX);

uPPR = @(x) kronPolyEval(K, x);
V = @(x) 0.5 * kronPolyEval(v,x);
dVdx = @(x) 0.5 * kronPolyDerivEval(v, x);

for i = 1:nY
    for j = 1:nX
        x = [X(i, j); Y(i, j)];
        valueFunction(i, j) = V(x);
        HJBResidual(i, j) = dVdx(x)*F(x) - 1/2*dVdx(x)*G(x)*G(x).'*dVdx(x).';
        if valueFunction(i, j) < 0
            valueFunction(i, j) = NaN;
        end
    end
end

fig1 = figure;
contourf(X, Y, valueFunction, 16, 'w'); hold on;
xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16); xticks([-pi, 0, pi]); xticklabels({'-\pi', '0', '\pi'})
load(fullfile('utils', 'YlGnBuRescaled.mat')); colormap(flip(YlGnBuRescaled))
% clim([0 4e4])
drawnow

axis off
fprintf('Exporting figure to: \n     plots/example26_valueFun_d%i_polynomial%i.pdf\n', degree, degree-1)
exportgraphics(fig1, sprintf('plots/example26_valueFun_d%i_polynomial%i.pdf', degree, degree-1), 'ContentType', 'vector', 'BackgroundColor', 'none');
colorbar

fig2 = figure;
pcolor(X, Y, log10(abs(HJBResidual))); shading interp;
xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16); xticks([-pi, 0, pi]); xticklabels({'-\pi', '0', '\pi'})
load('utils\YlGnBuRescaled.mat'); colormap(flip(YlGnBuRescaled))
clim([-3 9])
drawnow

axis off
fprintf('Exporting figure to: \n     plots/example26_valueFun-HJB-Error_d%i_polynomial%i.pdf\n', degree, degree-1)
exportgraphics(fig2, sprintf('plots/example26_valueFun-HJB-Error_d%i_polynomial%i.pdf', degree, degree-1), 'ContentType', 'vector', 'BackgroundColor', 'none');


%% Closed-loop phase portraits
plotPhasePortrait(degree, uPPR)

end



function plotPhasePortrait(degree, u)

% Write dynamics to file
fileID = fopen('examples\example26ExportedControlFunction.py', 'w');

% Define control law function declaration
funcStr = 'def u(z1, z2):\n    ';

% Now add the return statement where you print the actual control law in z
syms z1 z2
funcStr = [funcStr 'return ' char(vpa(u([z1;z2])))];

% Replace all the matlab stuff with python stuff
funcStr = strrep(funcStr, '^', '**'); % Replace MATLAB's ^ operator with Python's ** operator
funcStr = strrep(funcStr, 'sin', 'np.sin'); % Replace MATLAB's sin with np.sin
funcStr = strrep(funcStr, ';', ','); % Replace MATLAB's ; with , for the array

fprintf(fileID, funcStr);

fclose(fileID);

pyenv('ExecutionMode', 'OutOfProcess'); % Optional: You can switch to 'InProcess' mode if needed
pyrunfile("examples\plotExample26.py", controlDegree=degree-1)
open('plots/example26_phasePortrait.pdf')
terminate(pyenv);
end

