function [v, w] = runExample10(exportPlotData, varargin)
%runExample10 Runs 1D ODE example
%   Usage:  [v, w] = runExample10(exportPlotData)
%
%   Inputs:
%       exportPlotData     - Boolean variable to determine if plots/data are exported
%
%   Outputs:
%       v,w             are coefficients of the past and future energy
%                       function approximations, respectively.
%
%   References: [1] J. Borggaard and L. Zietsman, “Computation of nonlinear feedback
%   for flow control problems,” in 2018 American Control Conference (ACC), Jun. 2018,
%   doi: 10.23919/acc.2018.8431410.
%
%   Part of the NLbalancing repository.
%%

if nargin < 1
    exportPlotData = false;
end

xmax = 1;

[A, B, C, N, f, g, h] = getSystem10();
fprintf('Running Example 10\n')

eta = 1; % values should be between -\infty and 1.
% eta=1 is HJB/closed-loop balancing; this corresponds to reference [1]

fprintf('Simulating for eta=%g (gamma=%g)\n', eta, 1 / sqrt(1 - eta))

% Compute true energy functions
x = linspace(-xmax, xmax, 250);
EMinusAnalytic = EgammaMinusNumerical(x, f, g, h, eta);
EPlusAnalytic = EgammaPlusNumerical(x, f, g, h, eta);

%  Compute the polynomial approximations to the energy functions
d = 8;
[v] = approxPastEnergy(f, N, g, h, eta, d);
[w] = pqr(f, g, h, 1 / eta, d, true);

v2 = v{2}; v3 = v{3}; v4 = v{4}; v5 = v{5}; v6 = v{6}; v7 = v{7}; v8 = v{8};
w2 = w{2}; w3 = w{3}; w4 = w{4}; w5 = w{5}; w6 = w{6}; w7 = w{7}; w8 = w{8};

Ef2 = 0.5 * w2 * x .^ 2;
Ef3 = Ef2 + 0.5 * w3 * x .^ 3;
Ef4 = Ef3 + 0.5 * w4 * x .^ 4;
Ef5 = Ef4 + 0.5 * w5 * x .^ 5;
Ef6 = Ef5 + 0.5 * w6 * x .^ 6;
Ef7 = Ef6 + 0.5 * w7 * x .^ 7;
Ef8 = Ef7 + 0.5 * w8 * x .^ 8;

Ep2 = 0.5 * v2 * x .^ 2;
Ep3 = Ep2 + 0.5 * v3 * x .^ 3;
Ep4 = Ep3 + 0.5 * v4 * x .^ 4;
Ep5 = Ep4 + 0.5 * v5 * x .^ 5;
Ep6 = Ep5 + 0.5 * v6 * x .^ 6;
Ep7 = Ep6 + 0.5 * v7 * x .^ 7;
Ep8 = Ep7 + 0.5 * v8 * x .^ 8;

figure
subplot(2, 1, 1)
plot(x(1:10:end), EMinusAnalytic(1:10:end), '+', ...
    x, Ep2, ...
    x, Ep4, ...
    x, Ep6, ...
    x, Ep8, ...
    'LineWidth', 2)
legend('analytic', ...
    'degree 2', ...
    'degree 4', ...
    'degree 6', ...
    'degree 8', ...
    'Location', 'northeast')
xlabel('$x$', ...
    'interpreter', 'latex', ...
    'FontSize', 20, ...
    'fontweight', 'bold')
ylabel('$\mathcal{E}_\gamma^-$', ...
    'interpreter', 'latex', ...
    'FontSize', 20, ...
    'fontweight', 'bold')

xlim([-xmax xmax])
ylim([- .0025 .05])

subplot(2, 1, 2)
plot(x(1:10:end), EPlusAnalytic(1:10:end), '+', ...
    x, Ef2, ...
    x, Ef4, ...
    x, Ef6, ...
    x, Ef8, ...
    'LineWidth', 2)
legend('analytic', ...
    'degree 2', ...
    'degree 4', ...
    'degree 6', ...
    'degree 8', ...
    'Location', 'northeast')
xlabel('$x$', ...
    'interpreter', 'latex', ...
    'FontSize', 20, ...
    'fontweight', 'bold')
ylabel('$\mathcal{E}_\gamma^+$', ...
    'interpreter', 'latex', ...
    'FontSize', 20, ...
    'fontweight', 'bold')

xlim([-xmax xmax])
ylim([- .01 .1])

end

function [Ex] = EgammaPlusNumerical(xd, f, g, h, eta)

syms x;

a = -eta / 2 * (g{1}) ^ 2;
b = kronPolyEval(f, x);
c = 1/2 * kronPolyEval(h, x) ^ 2;

dEx2 = (-b - sqrt(b ^ 2 - 4 * a * c)) / (2 * a);
dEx1 = (-b + sqrt(b ^ 2 - 4 * a * c)) / (2 * a);

i = 1; Ex = zeros(1, length(xd));
for xi = xd
    if xi < 0
        Ex(i) = vpaintegral(dEx1, x, 0, xi);
    else
        Ex(i) = vpaintegral(dEx2, x, 0, xi);
    end
    i = i + 1;
end

end

function [Ex] = EgammaMinusNumerical(xd, f, g, h, eta)

syms x;

a = -1/2 * (g{1}) ^ 2;
b = -kronPolyEval(f, x);
c = eta / 2 * kronPolyEval(h, x) ^ 2;

dEx2 = (-b - sqrt(b ^ 2 - 4 * a * c)) / (2 * a);
dEx1 = (-b + sqrt(b ^ 2 - 4 * a * c)) / (2 * a);

i = 1; Ex = zeros(1, length(xd));
for xi = xd
    if xi < 0
        Ex(i) = vpaintegral(dEx1, x, 0, xi);
    else
        Ex(i) = vpaintegral(dEx2, x, 0, xi);
    end
    i = i + 1;
end
end
