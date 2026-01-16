function runExample11(nFterms, degree, reducedOrder)
%runExample11 Runs the 2D pendulum example. Plots a) the value function,
%   b) the HJB residual, and c) closed-loop phase portraits.
%
%   Usage:  runExample11(nFterms, degree)
%
%   runExample11() runs the default case of a degree 6 approximation.
%
%   Inputs:
%       nFterms         is the degree of dynamics
%       degree          is the degree of value function approximations
%       reducedOrder    reduced dimension for the higher-order coefficients
%
%   Generally one should use nFterms=degree-1, so that the value function
%   solution is the true Taylor expansion of the value function. If nFterms
%   is chosen to be smaller, local accuracy is sacrificed but global
%   behavior may be better due to avoiding "overfitting".
%
%   References:
%
%   Part of the PPR repository.
%% Process inputs
if nargin < 3
    if nargin < 2
        if nargin < 1
            nFterms = 5;
        end
        degree = nFterms + 1;
    end
    reducedOrder = 2;
end

fprintf('Running Example 11\n')

%% Get model
m = 1; L = 10; gravity = 9.81;
[f, g, ~, FofXU] = getSystem11(nFterms, m, L);

%% Open-loop phase portrait
% plotPhasePortrait(1, {[0, 0]})

%% Value function & HJB Residual plots
fprintf('Simulating...')
%  Compute the polynomial approximations to the past future energy function
% q = {0,0,sparse(linspace(1,2^3,2),1,0),sparse(linspace(1,2^4,2),1,100)}; r = {1,zeros(1,2)+0.0};
q = 0; r = 1;
options.verbose = true; options.reducedDimension = reducedOrder; tic
[v, K, ~] = ppr(f, g, q, r, degree, options);
fprintf("completed ppr() in %2.2f seconds. \n", toc)

nX = 301; nY = nX; xLim = pi; yLim = 5;
xPlot = linspace(-xLim, xLim, nX); yPlot = linspace(-yLim, yLim, nY); [X, Y] = meshgrid(xPlot, yPlot);

HJBResidual = zeros(nY, nX); valueFunction = zeros(nY, nX); 

for i = 1:nY
    for j = 1:nX
        x = [X(i, j); Y(i, j)];
        valueFunction(i, j) = 0.5 * kronPolyEval(v, x, degree=degree);
        HJBResidual(i, j) = computeHJBResidual(gravity, L, g, v, degree, x);
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
clim([0 4e4])
drawnow

% axis off
% fprintf('Exporting figure to: \n     plots/example11_valueFun_d%i_polynomial%i.pdf\n', degree, nFterms)
% exportgraphics(fig1, sprintf('plots/example11_valueFun_d%i_polynomial%i.pdf', degree, nFterms), 'ContentType', 'vector', 'BackgroundColor', 'none');
% colorbar

fig2 = figure;
pcolor(X, Y, log10(abs(HJBResidual))); shading interp;
xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16); xticks([-pi, 0, pi]); xticklabels({'-\pi', '0', '\pi'})
load('utils\YlGnBuRescaled.mat'); colormap(flip(YlGnBuRescaled))
clim([-3 9])
drawnow

% axis off
% fprintf('Exporting figure to: \n     plots/example11_valueFun-HJB-Error_d%i_polynomial%i.pdf\n', degree, nFterms)
% exportgraphics(fig2, sprintf('plots/example11_valueFun-HJB-Error_d%i_polynomial%i.pdf', degree, nFterms), 'ContentType', 'vector', 'BackgroundColor', 'none');


%% Closed-loop phase portraits
plotPhasePortrait(degree, K)

%% Evaluate region of attraction 
% In this section I would like to do some estimates of the region of
% attraction; there are a few approaches: 
%   1) Using the quadratic approximation of the closed-loop Lyapunov
%   function, can I prove a region, even if very conservative
%   2) Using the full approximation of the Lyapunov function, similar to
%   how I get the residual, can I get the region where the Lyapunov
%   function properties are satisfied? *** Note, the region D in Khalil Thm
%   4.1 is not the region of attraction! The largest lyapunov sublevel set
%   IN D is the estimate. 

% Method 2

end


% Helper functions
function c = estimateRegionOfAttraction(f,g,K,FofXU)
A_cl = f{1}-g{1}*K{1}; n = length(A_cl);

Q = eye(n);
P = lyap(A_cl.',Q); % P is positive definite for any positive definite Q


end



function plotPhasePortrait(degree, GainsPPR)
uPPR = @(z) (kronPolyEval(GainsPPR, z));

% Convert the function to a string
syms x1 x2
funcStr = char(vpa(uPPR([x1;x2]))); funcStr = strrep(funcStr, '^', '**'); % Replace MATLAB's ^ operator with Python's ** operator

% Create the Python function format
pythonFunc = ['def u(x1, x2):\n    return ' funcStr];

fileID = fopen('examples\example11ExportedControlFunction.py', 'w'); fprintf(fileID, pythonFunc); fclose(fileID);

pyenv('ExecutionMode', 'OutOfProcess'); % Optional: You can switch to 'InProcess' mode if needed
pyrunfile("examples\plotExample11.py", controlDegree=degree-1)
terminate(pyenv);
open('plots/example11_phasePortrait.pdf')
end

function res = computeHJBResidual(gravity, L, g, v, degree, x)

v = v(1:degree);

% constant input matrix B
res = (0.5 * kronPolyDerivEval(v, x)) * [x(2); 3 * gravity / (2 * L) * sin(x(1))] ...
    - 1 / 2 * 0.25 * kronPolyDerivEval(v, x) * g{1} * g{1}.' * kronPolyDerivEval(v, x).';

end
