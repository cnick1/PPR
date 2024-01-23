function runExample11_chisholmApproximant()
%runExample11_chisholmApproximant Runs the 2D example to plot energy functions as contour plots
%
%   Usage:  runExample11_chisholmApproximant()
%
%   References: [1]
%
%   Part of the NLbalancing repository.
exportPlotData = false; 
% close all; 

% Add Chisholm approximant matlab code to path
addpath('utils')

%% Get model and compute polynomial energy function
m = 1; L = 10; 
gravity = 9.81;
ell = 5; degree = ell+1; 
[f, g, h] = getSystem11(ell, m, L);
fprintf('Running Example 11: Chisholm approximant \n')

eta = 1;
[w] = pqr(f, g, h2q(h), eta, degree, true);

nX = 301; nY = nX;
xLim = pi; yLim = 5;
xPlot = linspace(-xLim, xLim, nX);
yPlot = linspace(-yLim, yLim, nY);
[X, Y] = meshgrid(xPlot, yPlot);

eFuture = zeros(nY, nX);
wRES = zeros(nY, nX);

for i = 1:nY
    for j = 1:nX
        x = [X(i, j); Y(i, j)];
        eFuture(i, j) = 0.5 * kronPolyEval(w, x, degree);
        wRES(i, j) = computeResidualFutureHJB_2D_example11(gravity, L, g, h, eta, w, degree, x);
        if eFuture(i, j) < 0
            eFuture(i, j) = NaN;
        end
    end
end

fig1 = figure;
contourf(X, Y, eFuture, 16, 'w'); hold on;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16)
xticks([-pi, 0, pi])
xticklabels({'-\pi', '0', '\pi'})
load(fullfile('utils', 'YlGnBuRescaled.mat'))
colormap(flip(YlGnBuRescaled))
clim([0 4e4])

fig2 = figure;
pcolor(X, Y, log10(abs(wRES))); shading interp;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16)
xticks([-pi, 0, pi])
xticklabels({'-\pi', '0', '\pi'})
colormap(flip(YlGnBuRescaled))
clim([-3 9])

if exportPlotData
    error("Need to update the file names")
    figure(fig1)
    axis off
    fprintf('Exporting figure to: \n     plots/example11_futureEnergy_d%i_polynomial%i.pdf\n', degree, nFterms)
    exportgraphics(fig1, sprintf('plots/example11_futureEnergy_d%i_polynomial%i.pdf', degree, nFterms), 'ContentType', 'vector', 'BackgroundColor', 'none');
    
    figure(fig2)
    axis off
    fprintf('Exporting figure to: \n     plots/example11_futureEnergy-HJB-Error_d%i_polynomial%i.pdf\n', degree, nFterms)
    exportgraphics(fig2, sprintf('plots/example11_futureEnergy-HJB-Error_d%i_polynomial%i.pdf', degree, nFterms), 'ContentType', 'vector', 'BackgroundColor', 'none');
end

figure(fig2)
colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex', 'XTick', -3:3:9, 'XTickLabel', {'1e-3', '1e0', '1e3', '1e6', '1e9'});
title(sprintf('Degree %i HJB Residual',degree))

figure(fig1)
colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex', 'XTick', 0:10000:4e4, 'XTickLabel', {'0', '1e4', '2e4', '3e4', '4e4'});
title(sprintf('Degree %i Future Energy Function',degree))

%% Compute a Chisholm approximant (rational approximant)
% Arrange energy function coefficients from w into C matrix for ACRS code
% to compute Chisholm approximant
C = zeros(15,15); % Just allocate plenty of space for C so that the code doesn't try to call and go outside of C 

n=2;

for k=[2,4]%,6]
    % Construct matrix idx where each row is the multi-index for one element of X
    idx = tt_ind2sub(ones(1, k) * n, (1:n ^ k)');

    for i = 1:n^k
        alpha = sum(idx(i,:) == 1);
        beta = sum(idx(i,:) == 2);
        C(alpha+1, beta+1) = C(alpha+1, beta+1) + w{k}(i);
    end
end


% Pick M and N and compute rational approximant w/ acrs fortran/matlab code
M = 3; N = 2; 
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
xticks([-pi, 0, pi])
xticklabels({'-\pi', '0', '\pi'})
colormap(flip(YlGnBuRescaled))
clim([0 4e4])

colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex', 'XTick', 0:10000:4e4, 'XTickLabel', {'0', '1e4', '2e4', '3e4', '4e4'});
title(sprintf('Rearranged Degree %i Future Energy Function; should reproduce fig1',degree))


% Plot Chisholm (rational) energy function given by R=A/B

Axy = evaluate2dPoly(A,X,Y);
Bxy = evaluate2dPoly(B,X,Y);
Rxy = 1/2.* Axy./Bxy;
Rxy(Rxy<0) = NaN;
% Rxy(Rxy>1e5) = NaN;

% figure
% % contourf(X,Y,Bxy,16)
% % surf(X,Y,Rxy,'edgecolor','none')
% contourf(X,Y,Rxy,16)
% colorbar

fig4 = figure;
contourf(X, Y, Rxy, levels, 'w'); hold on;
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16)
xticks([-pi, 0, pi])
xticklabels({'-\pi', '0', '\pi'})
colormap(flip(YlGnBuRescaled))
clim([0 4e4])

colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex', 'XTick', 0:10000:4e4, 'XTickLabel', {'0', '1e4', '2e4', '3e4', '4e4'});
title(sprintf('[%i,%i] Rational Future Energy Function',M,N))


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

function [res] = computeResidualFutureHJB_2D_example11(gravity, L, g, h, eta, w, degree, x)

w = w(1:degree);

%         constant B input
res = (0.5 * kronPolyDerivEval(w, x)) * [x(2); 3 * gravity / (2 * L) * sin(x(1))] ...
    - eta / 2 * 0.25 * kronPolyDerivEval(w, x) * g{1} * g{1}.' * kronPolyDerivEval(w, x).' ...
    + 0.5 * kronPolyEval(h, x).' * kronPolyEval(h, x);

%         Polynomial input
% res = (0.5 * kronPolyDerivEval(w, x)) * kronPolyEval(f, x) ...
%     - eta / 2 * 0.25 * kronPolyDerivEval(w, x) * (g{1} + kronPolyEval(g(2:end), x)) * (g{1} + kronPolyEval(g(2:end), x)).' * kronPolyDerivEval(w, x).' ...
%     + 0.5 * kronPolyEval(h, x).' * kronPolyEval(h, x);

end

function [N] = equivalenceClassIndices(n, k)
%equivalenceClassIndices - For a k-way tensor of dimension n (n^k entries), compute the matrix N which combines the equivalence class entries. This matrix essentially goes from Kronecker product form to the unique monomial form.
%
% Usage: [N] = equivalenceClassIndices(n,k)
%


%% Compute equivalence class index sets
% We will compute a vector like [1 5 5 7 5 7 7 8]; the number in the vector
% corresponds to the reference element for the equivalence class, and all
% of the entries with the same number are in the same equivalence class.
% For an n-dimensional k-order tensor (n^k entries), there are
% nchoosek(n+k-1,k) unique entries (distinct equivalence classes, i.e.
% monomials)

% Construct matrix ind where each row is the multi-index for one element of X
idx = tt_ind2sub(ones(1, k) * n, (1:n ^ k)');

% Find reference index for every element in the tensor - this is to its
% index in the symmetrized tensor. This puts every element into a 'class'
% of entries that will be the same under symmetry.
classidx = sort(idx, 2); % Normalize to one permutation, i.e. reference element
mult = [1 cumprod(ones(1, k - 1) * n)]; % Form shifts
linclassidx = (classidx - 1) * mult' + 1; % Form vector that maps to the reference elements

%% Form input-normal equations
% Get unique values from the input vector

diagIdxs = linspace(1, n ^ k, n); offdiagIdxs = setdiff(1:n ^ k, diagIdxs);
linclassidx(diagIdxs) = [];
[~, ~, uidx] = unique(linclassidx, 'stable');
Nhat = sparse(uidx, offdiagIdxs, 1, nchoosek(n + k - 1, k) - n, n ^ k);
N = [sparse(1:n, linspace(1, n ^ k, n), 1); Nhat];

end


