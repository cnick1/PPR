function runExample14_chisholmApproximant(ell,M,N)
%runExample14_chisholmApproximant Runs the 2D example to plot energy functions as contour plots
%
%   Usage:  runExample14_chisholmApproximant()
%
%   References: [1]
%
%   Part of the NLbalancing repository.
exportPlotData = false;
% close all;
if nargin < 3
    ell = 3; M = 3; N = 2;
end

% Add Chisholm approximant matlab code to path
addpath('utils')

%% Get model and compute polynomial energy function
% gravity = 9.81;
degree = ell+1;
[f, g, h] = getSystem14(ell,1);
fprintf('Running Example 14: Chisholm approximant \n')

% Compute open-loop OBSV energy
w = pqr(f, g, h2q(h), 0, degree, true);

nX = 301; nY = nX;
xLim = 3; yLim = 3;
xPlot = linspace(-xLim, xLim, nX);
yPlot = linspace(-yLim, yLim, nY);
[X, Y] = meshgrid(xPlot, yPlot);

eFuture = zeros(nY, nX);
% wRES = zeros(nY, nX);

for i = 1:nY
    for j = 1:nX
        x = [X(i, j); Y(i, j)];
        eFuture(i, j) = 0.5 * kronPolyEval(w, x, degree);
        % wRES(i, j) = computeResidualFutureHJB_2D_example14(gravity, L, g, h2q(h), 0, w, degree, x);
        if eFuture(i, j) < 0
            eFuture(i, j) = NaN;
        end
    end
end

load(fullfile('utils', 'YlGnBuRescaled.mat'))
if false
    fig1 = figure;
    contourf(X, Y, eFuture, 16, 'w'); hold on;
    xlabel('$x_1$', 'interpreter', 'latex');
    ylabel('$x_2$', 'interpreter', 'latex');
    set(gca, 'FontSize', 16)
    % xticks([-1, 0, 1])
    xticklabels({'1', '0', '1'})
    colormap(flip(YlGnBuRescaled))
    % clim([0 1.5])
    
    % fig2 = figure;
    % pcolor(X, Y, log10(abs(wRES))); shading interp;
    % xlabel('$x_1$', 'interpreter', 'latex');
    % ylabel('$x_2$', 'interpreter', 'latex');
    % set(gca, 'FontSize', 16)
    xticks([-1, 0, 1])
    xticklabels({'-1', '0', '1'})
    % colormap(flip(YlGnBuRescaled))
    % clim([-3 9])
    %
    % figure(fig2)
    % colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex', 'XTick', -3:3:9, 'XTickLabel', {'1e-3', '1e0', '1e3', '1e6', '1e9'});
    % title(sprintf('Degree %i HJB Residual',degree))
    
    figure(fig1)
    % colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex', 'XTick', 0:10000:4e4, 'XTickLabel', {'0', '1e4', '2e4', '3e4', '4e4'});
    title(sprintf('Degree %i Future Energy Function',degree))
end

%% Compute a Chisholm approximant (rational approximant)
% Arrange energy function coefficients from w into C matrix for ACRS code
% to compute Chisholm approximant
C = zeros(15,15); % Just allocate plenty of space for C so that the code doesn't try to call and go outside of C

n=2;

for k=2:degree
    % Construct matrix idx where each row is the multi-index for one element of X
    idx = tt_ind2sub(ones(1, k) * n, (1:n ^ k)');
    
    for i = 1:n^k
        alpha = sum(idx(i,:) == 1);
        beta = sum(idx(i,:) == 2);
        C(alpha+1, beta+1) = C(alpha+1, beta+1) + w{k}(i);
    end
end


% Pick M and N and compute rational approximant w/ acrs fortran/matlab code
% M = 3; N = 2;
[A, B] = acrs_sod(C, M, N);

% Display results
format shorte;

disp('Original C-matrix:');
disp(C(1:7,1:7));

disp('Resulting A-matrix:');
disp(A);
disp('Resulting B-matrix:');
disp(B);

%% Plot some preliminary results

% Plot polynomial energy function given by C; should match w from above
Cxy = 1/2 .* evaluate2dPoly(C,X,Y) ;% Energy function is 1/2*()
Cxy(Cxy<0) = NaN;

fig3 = figure;
[~, levels] = contourf(X, Y, Cxy, 16, 'w'); hold on;
levels = levels.LevelList;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16)
% xticks([-1, 0, 1])
% xticklabels({'-1', '0', '1'})
colormap(flip(YlGnBuRescaled))
% clim([0 1.5])

% colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex', 'XTick', 0:10000:4e4, 'XTickLabel', {'0', '1e4', '2e4', '3e4', '4e4'});
title(sprintf('Rearranged Degree %i Future Energy Function; should reproduce fig1',degree))


% Plot Chisholm (rational) energy function given by R=A/B

Axy = evaluate2dPoly(A,X,Y);
Bxy = evaluate2dPoly(B,X,Y);
Rxy = 1/2.* Axy./Bxy;
% Rxy(Rxy<0) = NaN;
% Rxy(Rxy>1e5) = NaN;


fig4 = figure;
contourf(X, Y, Rxy, levels, 'w'); hold on;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16)
% xticks([-1, 0, 1])
% xticklabels({'-1', '0', '1'})
colormap(flip(YlGnBuRescaled))
% clim([0 1.5])

if exportPlotData
    figure(fig4)
    axis off
    fprintf('Exporting figure to: \n     plots/example14_Chisholm_futureEnergy_m%i_n%i.pdf\n', M, N)
    exportgraphics(fig4, sprintf('plots/example14_Chisholm_futureEnergy_m%i_n%i.pdf', M, N), 'ContentType', 'vector', 'BackgroundColor', 'none');
end

% colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex', 'XTick', 0:10000:4e4, 'XTickLabel', {'0', '1e4', '2e4', '3e4', '4e4'});
title(sprintf('ℓ=%i, [%i,%i] Rational Future Energy Function',ell,M,N))
drawnow

%% Compute feedback laws
% Standard polynomial feedback
syms x1 x2
vpa((-g{1}.' * 0.5 * kronPolyDerivEval(w,[x1 ;x2]).'),2)

% Now compute it according to C(i,j)
C_dy = repmat(1:size(C,2)-1,size(C,2),1).*C(:,2:end);
U = -g{1}(2)*.5*C_dy(1:6,1:6);

u=evaluate2dPoly(U,x1,x2);

vpa(u,2)

% Rational feedback: R' = BA'-AB'/B^2
% Compute A' derivative wrt x2
A_dx2 = repmat(1:size(A,2)-1,size(A,2),1).*A(:,2:end);

% Compute B' derivative wrt x2
B_dx2 = repmat(1:size(B,2)-1,size(B,2),1).*B(:,2:end);

a=evaluate2dPoly(A,x1,x2);
b=evaluate2dPoly(B,x1,x2);
da=evaluate2dPoly(A_dx2,x1,x2);
db=evaluate2dPoly(B_dx2,x1,x2);

% ur = -g{1}(2)*.5* (b*da-a*db)/(b^2);
% vpa(expand(-g{1}(2)*.5*(b*da-a*db)),2)
% vpa(expand(b^2),2)

vpa(expand(-g{1}(2)*.5*(b*da-a*db))/expand(b^2))



end

function Pxy = evaluate2dPoly(C,X,Y)
% Given the coefficients of a polynomial f(x,y) = Sum_i,j C(i,j) x^i-1 y^j-1,
% and a mesh grid given by X,Y, return the values of f(x,y) on the grid.

Pxy = zeros(size(X));

for i=1:size(C,1)
    for j=1:size(C,2)
        Pxy = Pxy +  C(i,j) * X.^(i-1) .* Y.^(j-1);
    end
end

end

% function [res] = computeResidualFutureHJB_2D_example14(gravity, L, g, h, eta, w, degree, x)
%
% w = w(1:degree);
%
% %         constant B input
% res = (0.5 * kronPolyDerivEval(w, x)) * [x(2); 3 * gravity / (2 * L) * sin(x(1))] ...
%     - eta / 2 * 0.25 * kronPolyDerivEval(w, x) * g{1} * g{1}.' * kronPolyDerivEval(w, x).' ...
%     + 0;
% % + 0.5 * kronPolyEval(h, x).' * kronPolyEval(h, x);
%
%
% end
