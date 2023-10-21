%  This script runs two tests on the Burgers equation system.
%
%  The first is a computational performance benchmark with increasing system
%  order and fixed approximation degree.
%
%  The second uses a fixed system order and tests the accuracy of the
%  energy function approximations as the degree increases.
%
%  We assume this is run from the main directory and paths are set there.
%
%  This script generates the tables found in Section IV.C. of the paper
%    Nonlinear balanced truncation:  Part 1--Computing energy functions,
%    by Kramer, Gugercin, and Borggaard, arXiv:2209.07645.
%
%   Reference: [1] B. Kramer, S. Gugercin, J. Borggaard, and L. Balicki, “Nonlinear
%               balanced truncation: Part 1—computing energy functions,” arXiv,
%               Dec. 2022. doi: 10.48550/ARXIV.2209.07645
%              [2] B. Kramer, S. Gugercin, and J. Borggaard, “Nonlinear balanced
%               truncation: Part 2—model reduction on manifolds,” arXiv, Feb. 2023.
%               doi: 10.48550/ARXIV.2302.02036
%              [3] J. Borggaard and L. Zietsman, “On approximating polynomial-
%               -quadratic regulator problems,” IFAC-PapersOnLine, vol. 54, no. 9,
%               pp. 329–334, 2021, doi: 10.1016/j.ifacol.2021.06.090
%
%  Part of the NLbalancing and QQR repositories.
%%

eta = 0.9;

z_factor = 0.008; % scale factor on the initial conditions.
% zInit = 0.5*sin(2*pi*x)^2 on (0,0.5) and 0 otherwise
% so z_factor = 0.001 matches the initial condition in IV.C.

m = 4;
p = 4;

epsilon = 0.001;
alpha = 0.0; % alpha = 1.0, not controllable at n=8

%%
%  Computational performance of the energy function approximations.
%  Since the initial times are so short, we average nTest times
%
%  This builds TABLE II
%
fprintf('Table II Data\n')
nTest = 10;
degree = 3;

for n = [8, 16, 32, 64]
    fprintf('%d & ', n)
    [f, g, h, zInit] = getSystem3(n, m, p, epsilon, alpha);
    zInit = z_factor * zInit;

    tic; for i = 1:nTest, [w] = approxFutureEnergy(f, g, h, eta, degree); end, tt = toc / nTest;
    fprintf('%10.4e & ', length(w{degree}))
    fprintf('%8.2e & ', tt)

    for d = 2:degree, w{d} = w{d}.'; end
    wzInit = 0.5 * kronPolyEval(w, zInit, degree);
    fprintf('%12.6e \\\\ \n', wzInit)
end

%  Uncomment this block for the full table
nTest = 1;
degree = 3;

for n = [128, 256, 512, 1024]
    fprintf('%d & ', n)

    [f, g, h, zInit] = getSystem3(n, m, p, epsilon, alpha);
    zInit = z_factor * zInit;

    tic; for i = 1:nTest, [w] = approxFutureEnergy(f, g, h, eta, degree); end, tt = toc / nTest;
    fprintf('%10.4e & ', length(w{degree}))
    fprintf('%8.2e & ', tt)

    for d = 2:degree, w{d} = w{d}.'; end
    wzInit = 0.5 * kronPolyEval(w, zInit, degree);
    fprintf('%12.6e \\\\ \n', wzInit)
end

fprintf('\n\n')

%%
%  Computational performance of the energy function approximations.
%  Since the initial times are so short, we average nTest times
%
%  This builds TABLE III
%
fprintf('Table III Data\n')
nTest = 10;
degree = 4;

for n = [8, 16, 32, 64, 128]
    fprintf('%d & ', n)
    [f, g, h, zInit] = getSystem3(n, m, p, epsilon, alpha);
    zInit = z_factor * zInit;

    tic; for i = 1:nTest, [w] = approxFutureEnergy(f, g, h, eta, degree); end, tt = toc / nTest;
    fprintf('%10.4e & ', length(w{degree}))
    fprintf('%8.2e & ', tt)

    for d = 2:degree, w{d} = w{d}.'; end
    wzInit = 0.5 * kronPolyEval(w, zInit, degree);
    fprintf('%12.6e \\\\ \n', wzInit)
end

fprintf('\n\n')

%%
%  Test convergence of the energy function at a "point" z_factor*zInit
%  with increasing degree...
%
%  This builds TABLE IV
%
fprintf('Table IV Data\n')
n = 8;
[f, g, h, zInit] = getSystem3(n, m, p, epsilon, alpha);
zInit = z_factor * zInit;
for degree = [2, 3, 4, 5, 6, 7, 8]
    fprintf('%d & ', degree)

    [v] = approxPastEnergy(f, g, h, eta, degree);
    for d = 2:degree, v{d} = v{d}.'; end
    vzInit = 0.5 * kronPolyEval(v, zInit, degree);
    fprintf('%12.6e & ', vzInit)

    [w] = approxFutureEnergy(f, g, h, eta, degree);
    for d = 2:degree, w{d} = w{d}.'; end
    wzInit = 0.5 * kronPolyEval(w, zInit, degree);
    fprintf('%12.6e \\\\ \n', wzInit)
end

%   %% We can check the eigenvalues of A as one test that matrices are
%   %  computed correctly.
%   g = sort(eig(A)','descend'); disp(g(1:8))
%
%   % should converge to the following as n increases:
%   disp(-epsilon*pi^2*(1:8).^2+alpha)
