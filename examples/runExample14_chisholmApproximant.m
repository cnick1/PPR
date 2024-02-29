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
load(fullfile('utils', 'YlGnBuRescaled.mat'))

%% Get model and compute polynomial energy function
degree = ell+1;
[f, g, h] = getSystem14(ell,1);
fprintf('Running Example 14: Chisholm approximant \n')

% Compute open-loop OBSV energy
q = h2q(h); R = 0;
options.verbose = true; options.skipGains = true;
w = ppr(f, g, q, R, degree, options);
W2 = reshape(w{2},2,2); [V, D] = eig(W2);

%% Compute a Chisholm approximant (rational approximant)
% Arrange energy function coefficients from w into C matrix for ACRS code
% to compute Chisholm approximant
C = kronCoeffs2coeffTable(w);

% Pick M and N and compute rational approximant w/ acrs fortran/matlab code
[A, B] = acrs_sod(C, M, N);
a_coeff = coeffTable2kronCoeffs(A); b_coeff = coeffTable2kronCoeffs(B);

% Display results
format shorte;

disp('Original C-matrix:'); disp(C);

disp('Resulting A-matrix:'); disp(A);
disp('Resulting B-matrix:'); disp(B);

%% Plot some preliminary results
xLim = 1; yLim = 1;
nX = 301; nY = nX;
xPlot = linspace(-xLim, xLim, nX);
yPlot = linspace(-yLim, yLim, nY);
[X, Y] = meshgrid(xPlot, yPlot);

poly_value = zeros(nY, nX); rat_value = zeros(nY, nX);
poly_RES = zeros(nY, nX); rat_RES = zeros(nY, nX);

% Compute value functions and HJB errors
for i = 1:nY
    for j = 1:nX
        x = [X(i, j); Y(i, j)];
        
        poly_value(i, j) = 0.5 * kronPolyEval(w, x);
        poly_RES(i, j) = computeResidualFutureHJB_2D_example14(q, R, 0.5*kronPolyDerivEval(w, x), x);
        
        a = kronPolyEval(a_coeff, x); b = 1+kronPolyEval(b_coeff, x);
        da = kronPolyDerivEval(a_coeff, x); db = kronPolyDerivEval(b_coeff, x);
        
        rat_value(i, j) = 0.5 * a/b;
        rat_RES(i, j) = computeResidualFutureHJB_2D_example14(q, R, 0.5*(b*da-a*db)/b^2, x);
    end
end

poly_value(poly_value<0) = NaN;
rat_value(rat_value<0) = NaN;

fig1 = figure;
[~, levels] = contourf(X, Y, poly_value, 16, 'w'); hold on;
levels = levels.LevelList;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16)
colormap(flip(YlGnBuRescaled))
% clim([0 4e4])

quiver(0,0,V(1,1),V(2,1))

colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex');
title(sprintf('Rearranged Degree %i Future Energy Function; should reproduce fig1',degree))

fig2 = figure;
contourf(X, Y, rat_value, levels, 'w'); hold on;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16)
colormap(flip(YlGnBuRescaled))
% clim([0 4e4])

if exportPlotData
    figure(fig2)
    axis off
    fprintf('Exporting figure to: \n     plots/example14_Chisholm_futureEnergy_m%i_n%i.pdf\n', M, N)
    exportgraphics(fig2, sprintf('plots/example14_Chisholm_futureEnergy_m%i_n%i.pdf', M, N), 'ContentType', 'vector', 'BackgroundColor', 'none');
end

quiver(0,0,V(1,1),V(2,1))

colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex');
title(sprintf('ℓ=%i, [%i,%i] Rational Future Energy Function',ell,M,N))
drawnow

fig3 = figure;
pcolor(X, Y, log10(abs(poly_RES))); shading interp; hold on;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16)
colormap(flip(YlGnBuRescaled))
% clim([-3 9])
quiver(0,0,V(1,1),V(2,1))

figure(fig3)
colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex');
title(sprintf('Degree %i HJB Residual',degree))

fig4 = figure;
pcolor(X, Y, log10(abs(rat_RES))); shading interp; hold on;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16)
colormap(flip(YlGnBuRescaled))
% clim([-3 9])
quiver(0,0,V(1,1),V(2,1))

figure(fig4)
colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex');
title(sprintf('R[%i,%i] HJB Residual',M,N))

%% Compute feedback laws
% Standard polynomial feedback
syms x1 x2

u_poly = -g{1}.'*0.5*kronPolyDerivEval(w,[x1 ;x2]).';

% Rational feedback: R' = BA'-AB'/B^2
a = kronPolyEval(a_coeff, [x1 ;x2]); b = 1+kronPolyEval(b_coeff, [x1 ;x2]);
da = kronPolyDerivEval(a_coeff, [x1 ;x2]).'; db = kronPolyDerivEval(b_coeff, [x1 ;x2]).';

u_rat = -g{1}.'*.5*expand(b*da-a*db)/expand(b^2);

vpa(u_poly,2)
vpa(u_rat,2)

end

function [res] = computeResidualFutureHJB_2D_example14(q, R, dVdx, x)
fx = [(10*sin(x(2) + x(1))*(cos(x(2)) + 1))/(2*cos(x(2)) - (cos(x(2)) + 1)^2 + 3) - (10*sin(x(2) + x(1)) + 20*sin(x(1)))/(2*cos(x(2)) - (cos(x(2)) + 1)^2 + 3);
    ((cos(x(2)) + 1)*(10*sin(x(2) + x(1)) + 20*sin(x(1))))/(2*cos(x(2)) - (cos(x(2)) + 1)^2 + 3) - (10*sin(x(2) + x(1))*(2*cos(x(2)) + 3))/(2*cos(x(2)) - (cos(x(2)) + 1)^2 + 3)];
gx = [1/(2*cos(x(2)) - (cos(x(2)) + 1)^2 + 3);
    -(cos(x(2)) + 1)/(2*cos(x(2)) - (cos(x(2)) + 1)^2 + 3)];
qx = kronPolyEval(q,x);

if R == 0
    res = dVdx * fx + 1/2*qx;
else
    res = dVdx * fx - 1/2 * dVdx * gx * (R \ gx.' * dVdx.') + 1/2*qx;
end
end

