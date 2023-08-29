function [v, w] = runExample5(exportPlotData, nGterms, degree, eps, eta, varargin)
%runExample5 Runs the 3D unicycle example
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
%
%       \dot(x1) = u1 cos x3
%       \dot(x2) = u1 sin x3  - eps*x2
%       \dot(x3) = u2
%             y = x
%
%   References: [1]
%
%   Part of the NLbalancing repository.
%% Process inputs
if nargin < 5
    if nargin < 4
        if nargin < 3
            if nargin < 2
                if nargin < 1
                    exportPlotData = 0;
                end
                nGterms = 7;
            end
            degree = 2;
        end
        eps = 0.1;
    end
    % Compute energy functions
    eta = 1; % values should be between -\infty and 1.
    % eta=1 is HJB/closed-loop balancing, 0 is open loop.
end

load(fullfile('utils', 'YlGnBuRescaled.mat'));

%% Get model and compute energy functions
[A, ~, C, N, f, g, h] = getSystem5(eps, nGterms);
fprintf('Running Example 5\n')
fprintf('Simulating for eta=%g (gamma=%g)\n', eta, 1 / sqrt(1 - eta))

[w] = approxFutureEnergy(f, N, g, h, eta, degree, true);

%% 1D plots
xrange = 1; nPoints = 100;

figure;

% Iterate through each axis and plot
n = length(A);
for xi = 1:n
    xs = zeros(nPoints,n);
    xs(:,xi) = linspace(-xrange,xrange,nPoints).';
    
    energies = zeros(nPoints,1);
    for i=1:nPoints
        x = xs(i,:).';
        energies(i) = 0.5*kronPolyEval(w, x, degree);
    end
    
    subplot(1,n,xi)
    plot(xs(:,xi), energies,'LineWidth', 2)
    
    xlabel(sprintf('$x_%i$',xi), 'interpreter', 'latex', 'FontSize', 20, 'fontweight', 'bold'); ylabel('$\mathcal{E}_\gamma^+$', 'interpreter', 'latex', 'FontSize', 20,  'fontweight', 'bold')
    
end


%% 2D plots
xPlot = linspace(-xrange, xrange, nPoints);
yPlot = linspace(-xrange, xrange, nPoints);
[X, Y] = meshgrid(xPlot, yPlot);

eFuture = zeros(nPoints, nPoints);
wRES = zeros(nPoints, nPoints);

insert = @(a, x, n)cat(2,  x(1:n-1).', a, x(n:end).').';

figure('Position',[63.6667 1 1644 887.3333]);
for xi = 1:n
    for i = 1:nPoints
        for j = 1:nPoints
            x = [X(i, j); Y(i, j)];
            x = insert(0,x,xi);
            eFuture(i, j) = 0.5 * kronPolyEval(w, x, degree);
            if eFuture(i, j) < 0
                eFuture(i, j) = NaN;
            end
            wRES(i,j) = computeResidualFutureHJB_2D_example5(f, g, h, eta, w, degree, x);
        end
    end
    
    subplot(2,n,n-xi+1);
    %     if xi == 1
    %         contourf(X, Y, eFuture, [0 1.7 3.5 5.3 7.0 8.7 10.4 12.2 14.0 15.7 17.4 19.2 20.9 22.6 24.4], 'w');
    %         caxis([0 28])
    %     elseif xi == 2
    %         contourf(X, Y, eFuture, [0 1.7 3.5 5.3 7.0 8.7 10.4 12.2 14.0 15.7 17.4 19.2 20.9 22.6 24.4]./3, 'w');
    % %         contourf(X, Y, eFuture, 16, 'w','ShowText','on');
    %         caxis([0 9])
    %     elseif xi == 3
    %         contourf(X, Y, eFuture, [0 1.7 3.5 5.3 7.0 8.7 10.4 12.2 14.0 15.7 17.4 19.2 20.9 22.6 24.4], 'w');
    %         caxis([0 28])
    %     end
    if xi == 1
        contourf(X, Y, eFuture, [0 1.7 3.5 5.3 7.0 8.7 10.4 12.2 14.0 15.7 17.4 19.2 20.9 22.6 24.4]./10, 'w');
        caxis([0 2.8])
        xlabel('$x_2$', 'interpreter', 'latex'); ylabel('$x_3$', 'interpreter', 'latex');
    elseif xi == 2
        contourf(X, Y, eFuture, [0 1.7 3.5 5.3 7.0 8.7 10.4 12.2 14.0 15.7 17.4 19.2 20.9 22.6 24.4]./30, 'w');
        %         contourf(X, Y, eFuture, 16, 'w','ShowText','on');
        caxis([0 .9])
        xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_3$', 'interpreter', 'latex');
        
    elseif xi == 3
        contourf(X, Y, eFuture, [0 1.7 3.5 5.3 7.0 8.7 10.4 12.2 14.0 15.7 17.4 19.2 20.9 22.6 24.4]./10, 'w');
        caxis([0 2.8])
        xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex');
    end
    colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex'); set(gca, 'FontSize', 16)
    %     xticks([-xrange,0, xrange]); % xticklabels({'-\pi','0','\pi'})
    axis equal; colormap(flip(YlGnBuRescaled))
    
    
    % HJB Residual Plots
    subplot(2,n,n-xi+n+1);
    %contourf(X, Y, abs(wRES), 16, 'w'); hold on;
    pcolor(X, Y, log10(abs(wRES))); shading interp; colorbar;
    xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex');
    set(gca, 'FontSize', 16)
    colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex','XTick', -2:4,'XTickLabel',{'1e-2','1e-1','1e0','1e1','1e2','1e3','1e4'});
    %     xticks([-xrange,0, xrange]); % xticklabels({'-\pi','0','\pi'})
    %     xlim([-xrange, xrange]);ylim([-xrange, xrange]);
    colormap(flip(YlGnBuRescaled))
    caxis([-2 4])
    
end




end


function [res] = computeResidualFutureHJB_2D_example5(f, g, h, eta, w, degree, x)

w = w(1:degree);
%       \dot(x1) = u1 cos x3
%       \dot(x2) = u1 sin x3  - x2
%       \dot(x3) = u2

%         Polynomial input
res = (0.5 * kronPolyDerivEval(w, x)) * kronPolyEval(f, x) ...
    - eta / 2 * 0.25 * kronPolyDerivEval(w, x) * ([cos(x(3)),0;sin(x(3)),0;0,1]) * ([cos(x(3)),0;sin(x(3)),0;0,1]).' * kronPolyDerivEval(w, x).' ...
    + 0.5 * kronPolyEval(h, x).' * kronPolyEval(h, x);

end

