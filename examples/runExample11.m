function [w] = runExample11(nFterms, degree, varargin)
%runExample11 Runs the 2D example to plot energy functions as contour plots
%
%   Usage:  [v,w] = runExample11(degree,plotEnergy,plotBalancing,balancingDegree,
%                               numGTermsModel, numGTermsApprox, exportPlotData, kawanoModel)
%
%   runExample11() runs the default case of a quadratic model from [1] which
%                 is based on a model from [2].
%
%   Inputs:
%       degree          is the degree of energy function approximations
%       exportPlotData   Boolean variable to determine if plots/data are exported
%
%
%   Outputs:
%       v,w              are coefficients of the past and future energy
%                        function approximations, respectively.
%
%   The value of eta is set below.
%
%   References: [1]
%
%   Part of the NLbalancing repository.
%% Process inputs
if nargin < 2
    if nargin < 1
        nFterms = 5;
    end
    degree = 8;
end

if nFterms == 1
    nFterms = 2; % Note F2 is zero, this is just to be able to compute a controller and ignore the error if F2 doesn't exist
end
fprintf('Running Example 11\n')

%% Get model
m = 1; L = 10;
gravity = 9.81;
[f, g, ~] = getSystem11(nFterms, m, L);

%% Open-loop phase portrait
% Create the Python function format
pythonFunc = ['def u(x1, x2):\n    return 0' ];
% Write the function to a Python file
fileID = fopen('examples\example11ExportedControlFunction.py', 'w'); fprintf(fileID, pythonFunc); fclose(fileID);

pyenv('ExecutionMode', 'OutOfProcess'); % Optional: You can switch to 'InProcess' mode if needed
pyrunfile("examples\plotExample11.py", controlDegree=0)
terminate(pyenv);
open('plots/example11_phasePortrait.pdf')

%% Energy function plots
fprintf('Simulating...')

%  Compute the polynomial approximations to the past future energy function
options.verbose = true; tic
[w, GainsPPR, ~] = ppr(f, g, 0, 1, degree, options);
fprintf("completed ppr() in %2.2f seconds. \n", toc)

nX = 301; nY = nX; xLim = pi; yLim = 5;
xPlot = linspace(-xLim, xLim, nX); yPlot = linspace(-yLim, yLim, nY); [X, Y] = meshgrid(xPlot, yPlot);

wRES = zeros(nY, nX); eFuture = zeros(nY, nX);

for i = 1:nY
    for j = 1:nX
        x = [X(i, j); Y(i, j)];
        eFuture(i, j) = 0.5 * kronPolyEval(w, x, degree);
        wRES(i, j) = computeResidualFutureHJB_2D_example11(gravity, L, g, w, degree, x);
        if eFuture(i, j) < 0
            eFuture(i, j) = NaN;
        end
    end
end

fig1 = figure;
contourf(X, Y, eFuture, 16, 'w'); hold on;
xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16); xticks([-pi, 0, pi]); xticklabels({'-\pi', '0', '\pi'})
load(fullfile('utils', 'YlGnBuRescaled.mat')); colormap(flip(YlGnBuRescaled))

fig2 = figure;
pcolor(X, Y, log10(abs(wRES))); shading interp;
xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16); xticks([-pi, 0, pi]); xticklabels({'-\pi', '0', '\pi'})
load('utils\YlGnBuRescaled.mat'); colormap(flip(YlGnBuRescaled))
caxis([-3 9])
drawnow 

%% Closed-loop phase portraits
uPPR = @(z) (kronPolyEval(GainsPPR, z));

% Convert the function to a string
syms x1 x2
funcStr = char(vpa(uPPR([x1;x2]))); funcStr = strrep(funcStr, '^', '**'); % Replace MATLAB's ^ operator with Python's ** operator

% Create the Python function format
pythonFunc = ['def u(x1, x2):\n    return ' funcStr];

% Write the function to a Python file
fileID = fopen('examples\example11ExportedControlFunction.py', 'w'); fprintf(fileID, pythonFunc); fclose(fileID);

pyenv('ExecutionMode', 'OutOfProcess'); % Optional: You can switch to 'InProcess' mode if needed
pyrunfile("examples\plotExample11.py", controlDegree=degree-1)
terminate(pyenv);
open('plots/example11_phasePortrait.pdf')
end

function [res] = computeResidualFutureHJB_2D_example11(gravity, L, g, w, degree, x)

w = w(1:degree);

%         constant B input
res = (0.5 * kronPolyDerivEval(w, x)) * [x(2); 3 * gravity / (2 * L) * sin(x(1))] ...
    - 1 / 2 * 0.25 * kronPolyDerivEval(w, x) * g{1} * g{1}.' * kronPolyDerivEval(w, x).';


end
