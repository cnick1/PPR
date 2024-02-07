function runExample12_chisholmApproximant(ell,M,N)
%runExample12_chisholmApproximant Runs the 2D example to plot energy functions as contour plots
%
%   Usage:  runExample12_chisholmApproximant()
%
%   References: [1]
%
%   Part of the NLbalancing repository.
exportPlotData = false;
% close all;
if nargin < 3
    ell = 5; M = 2; N = 1;
end

% Add Chisholm approximant matlab code to path
addpath('utils')
load(fullfile('utils', 'YlGnBuRescaled.mat'))

%% Get model and compute polynomial energy function
degree = ell+1;
[f, g, h] = getSystem12(ell, false);
fprintf('Running Example 12: Chisholm approximant \n')

% Compute open-loop OBSV energy
q = h2q(h); R = 0;
w = ppr(f, g, q, R, degree, true);
W2 = reshape(w{2},2,2); [V, D] = eig(W2);

%% Compute a Chisholm approximant (rational approximant)
% Arrange energy function coefficients from w into C matrix for ACRS code
% to compute Chisholm approximant
C = kronCoeffs2coeffTable(w);

format shorte;
disp('Original C-matrix:'); disp(C);

% Pick M and N and compute rational approximant w/ acrs fortran/matlab code
[A, B] = acrs_sod(C, M, N);
A(isnan(A)) = 0; B(isnan(B)) = 0;
a_coeff = coeffTable2kronCoeffs(A); b_coeff = coeffTable2kronCoeffs(B);

% Display results
disp('Resulting A-matrix:'); disp(A);
disp('Resulting B-matrix:'); disp(B);

%% Plot energy functions and HJB residual errors
xLim = 1; yLim = 1;
nX = 301; nY = nX;
xPlot = linspace(-xLim, xLim, nX);
yPlot = linspace(-yLim, yLim, nY);
[X, Y] = meshgrid(xPlot, yPlot);

poly_value = zeros(nY, nX); rat_value = zeros(nY, nX);
an_value = zeros(nY, nX); an_RES = zeros(nY, nX);
poly_RES = zeros(nY, nX); rat_RES = zeros(nY, nX);

% Compute value functions and HJB errors
for i = 1:nY
    for j = 1:nX
        x = [X(i, j); Y(i, j)];
        
        an_value(i, j) = 1/2*(36*x(1)^2 + 9*x(2)^2 + 18*x(1)^3*x(2) + 18*x(1)*x(2)^3 + x(1)^6 + 6*x(1)^4*x(2)^2 + 9*x(1)^2*x(2)^4 + 4*x(2)^6)/(1 + x(1)^4 + 2*x(1)^2*x(2)^2 + x(2)^4);
        
        EanalyticGradient = [(9*x(2)^3 + 36*x(1) + 9*x(1)*x(2)^4 + 27*x(1)^2*x(2) + 12*x(1)^3*x(2)^2 + 3*x(1)^5)/(1 + x(2)^4 + 2*x(1)^2*x(2)^2 + x(1)^4) - ((4*x(1)*x(2)^2 + 4*x(1)^3)*((9*x(2)^2)/2 + 2*x(2)^6 + 9*x(1)*x(2)^3 + 18*x(1)^2 + (9*x(1)^2*x(2)^4)/2 + 9*x(1)^3*x(2) + 3*x(1)^4*x(2)^2 + x(1)^6/2))/(2*x(1)^2*x(2)^2 + x(1)^4 + x(2)^4 + 1)^2;
            (9*x(2) + 12*x(2)^5 + 27*x(1)*x(2)^2 + 18*x(1)^2*x(2)^3 + 9*x(1)^3 + 6*x(1)^4*x(2))/(1 + x(2)^4 + 2*x(1)^2*x(2)^2 + x(1)^4) - ((4*x(2)^3 + 4*x(1)^2*x(2))*((9*x(2)^2)/2 + 2*x(2)^6 + 9*x(1)*x(2)^3 + 18*x(1)^2 + (9*x(1)^2*x(2)^4)/2 + 9*x(1)^3*x(2) + 3*x(1)^4*x(2)^2 + x(1)^6/2))/(2*x(1)^2*x(2)^2 + x(1)^4 + x(2)^4 + 1)^2];
        an_RES(i, j) = computeResidualFutureHJB_2D_example12(q, R, EanalyticGradient.', x);

        poly_value(i, j) = 0.5 * kronPolyEval(w, x);
        poly_RES(i, j) = computeResidualFutureHJB_2D_example12(q, R, 0.5 * kronPolyDerivEval(w, x), x);
        
        a = kronPolyEval(a_coeff, x); b = 1+kronPolyEval(b_coeff, x);
        da = kronPolyDerivEval(a_coeff, x); db = kronPolyDerivEval(b_coeff, x);
        
        rat_value(i, j) = 0.5 * a/b;
        rat_RES(i, j) = computeResidualFutureHJB_2D_example12(q, R, 0.5*(b*da-a*db)/b^2, x);
    end
end

poly_value(poly_value<0) = NaN;
rat_value(rat_value<0) = NaN;

% Plot analytical energy function
fig1p5 = figure;
[~, levels] = contourf(X, Y, an_value, 16, 'w'); hold on;
levels = levels.LevelList;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16)
colormap(flip(YlGnBuRescaled))
% clim([0 4e4])

quiver(0,0,V(1,1),V(2,1))

colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex');
title('Analytical Future Energy Function')

fig1 = figure;
contourf(X, Y, poly_value, levels, 'w'); hold on;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16)
colormap(flip(YlGnBuRescaled))
% clim([0 4e4])

quiver(0,0,V(1,1),V(2,1))

colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex');
title(sprintf('Degree %i Future Energy Function',degree))

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
    fprintf('Exporting figure to: \n     plots/example12_Chisholm_futureEnergy_m%i_n%i.pdf\n', M, N)
    exportgraphics(fig2, sprintf('plots/example12_Chisholm_futureEnergy_m%i_n%i.pdf', M, N), 'ContentType', 'vector', 'BackgroundColor', 'none');
end

quiver(0,0,V(1,1),V(2,1))

colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex');
title(sprintf('ℓ=%i, [%i,%i] Rational Future Energy Function',ell,M,N))
drawnow

fig3p5 = figure;
pcolor(X, Y, log10(abs(an_RES))); shading interp; hold on;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16)
colormap(flip(YlGnBuRescaled))
clim([-10 2])
quiver(0,0,V(1,1),V(2,1))

figure(fig3p5)
colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex');
title(sprintf('Analytical HJB Residual (should be zero)'))

fig3 = figure;
pcolor(X, Y, log10(abs(poly_RES))); shading interp; hold on;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16)
colormap(flip(YlGnBuRescaled))
clim([-10 2])
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
clim([-10 2])
quiver(0,0,V(1,1),V(2,1))

figure(fig4)
colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex');
title(sprintf('R[%i,%i] HJB Residual',M,N))


figure(fig2)



%% Compute feedback laws
% Standard polynomial feedback
% syms x1 x2
% 
% u_poly = -g{1}.'*0.5*kronPolyDerivEval(w,[x1 ;x2]).';
% 
% % Rational feedback: R' = BA'-AB'/B^2
% a = kronPolyEval(a_coeff, [x1 ;x2]); b = 1+kronPolyEval(b_coeff, [x1 ;x2]);
% da = kronPolyDerivEval(a_coeff, [x1 ;x2]).'; db = kronPolyDerivEval(b_coeff, [x1 ;x2]).';
% 
% u_rat = -g{1}.'*.5*expand(b*da-a*db)/expand(b^2);
% 
% vpa(u_poly,2)
% vpa(u_rat,2)

end

function [res] = computeResidualFutureHJB_2D_example12(q, R, dVdx, x)
fx = [-9 * x(1) + 6 * x(1) ^ 2 * x(2) + 6 * x(2) ^ 3 - x(1) ^ 5 - 2 * x(1) ^ 3 * x(2) ^ 2 - x(1) * x(2) ^ 4;
    -9 * x(2) - 6 * x(1) ^ 3 - 6 * x(1) * x(2) ^ 2 - x(1) ^ 4 * x(2) - 2 * x(1) ^ 2 * x(2) ^ 3 - x(2) ^ 5];
gx = [3 * sqrt(2) * (9 - 6 * x(1) * x(2) + x(1) ^ 4 - x(2) ^ 4) / (9 + x(1) ^ 4 + 2 * x(1) ^ 2 * x(2) ^ 2 + x(2) ^ 4), sqrt(2) * (-9 * x(1) ^ 2 - 27 * x(2) ^ 2 + 6 * x(1) ^ 3 * x(2) + 6 * x(1) * x(2) ^ 3 - (x(1) ^ 2 + x(2) ^ 2) ^ 3) / (9 + x(1) ^ 4 + 2 * x(1) ^ 2 * x(2) ^ 2 + x(2) ^ 4);
    sqrt(2) * (27 * x(1) ^ 2 + 9 * x(2) ^ 2 + 6 * x(1) ^ 3 * x(2) + 6 * x(1) * x(2) ^ 3 + (x(1) ^ 2 + x(2) ^ 2) ^ 3) / (9 + x(1) ^ 4 + 2 * x(1) ^ 2 * x(2) ^ 2 + x(2) ^ 4), 3 * sqrt(2) * (9 + 6 * x(1) * x(2) - x(1) ^ 4 + x(2) ^ 4) / (9 + x(1) ^ 4 + 2 * x(1) ^ 2 * x(2) ^ 2 + x(2) ^ 4)];
% qx = kronPolyEval(q,x);
hx = [2 * sqrt(2) * (3 * x(1) + x(1) ^ 2 * x(2) + x(2) ^ 3) * (3 - x(1) ^ 4 - 2 * x(1) ^ 2 * x(2) ^ 2 - x(2) ^ 4) / (1 + x(1) ^ 4 + 2 * x(1) ^ 2 * x(2) ^ 2 + x(2) ^ 4);
            sqrt(2) * (3 * x(2) - x(1) ^ 3 - x(1) * x(2) ^ 2) * (3 - x(1) ^ 4 - 2 * x(1) ^ 2 * x(2) ^ 2 - x(2) ^ 4) / (1 + x(1) ^ 4 + 2 * x(1) ^ 2 * x(2) ^ 2 + x(2) ^ 4)];
qx = hx.'*hx;

if R == 0
    res = dVdx * fx + 1/2*qx;
else
    res = dVdx * fx - 1/2*R * dVdx * gx * (R \ gx.' * dVdx.') + 1/2*qx;
end
end