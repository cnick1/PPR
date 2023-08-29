function [v, w] = runExample7_reproduceResults()
%runExample7 Runs 3D aircraft stall model
%   Usage:  [v, w] = runExample7()
%
%   runExample7() runs TODO.
%
%   Outputs:
%       v,w             are coefficients of the past and future energy
%                       function approximations, respectively.
%
%   Part of the NLbalancing repository.
%%
[~, ~, ~, ~, f, g, h] = getSystem7();
% g = g(1);
fprintf('Running Example 7\n')

%% Define original dynamics and controllers
F = @(x) kronPolyEval(f, x);

% This is the reduced dynamics
g = g{1};
G = @(x) (g);
U2 = @(x) [0; 0; 0]; U3 = [0; 0; 0];

% This is technically the full model that Garrard presented
% G = @(x) (g{1} + kronPolyEval(g(2:end), x));
% U2 = @(x) [.47*x(1);0;46]; U3 = [.63;0;61.4];

%% Find stall angle
tspan = [0, 10];

% First IC
alpha0 = 31.6;
X0 = [pi / 179 * alpha0; 0; 0];

[t1, X1] = ode45(@(t, x) F(x), tspan, X0);

figure; hold on;
plot(t1, X1(:, 1) / (pi / 180))
plot(t1, X1(:, 2) / (pi / 180))
plot(t1, X1(:, 3) / (pi / 180))

%% Recreate Garrard Figure 1
tspan = [0, 5];

uLinGarrard = @(x) (- 0.053 * x(1) + 0.5 * x(2) + 0.521 * x(3));
uQuadGarrard = @(x) (- 0.053 * x(1) + 0.5 * x(2) + 0.521 * x(3) + 0.04 * x(1) ^ 2 - 0.048 * x(1) * x(2));
uCubGarrard = @(x) (- 0.053 * x(1) + 0.5 * x(2) + 0.521 * x(3) + 0.04 * x(1) ^ 2 - 0.048 * x(1) * x(2) + 0.374 * x(1) ^ 3 - 0.312 * x(1) ^ 2 * x(2));

% First IC
alpha0 = 22.9;
X0 = [pi / 180 * alpha0; 0; 0];

u = uLinGarrard; [t1, X1] = ode45(@(t, x) F(x) + G(x) * u(x) + U2(x) * u(x) ^ 2 + U3 * u(x) ^ 3, tspan, X0);
u = uQuadGarrard; [t2, X2] = ode45(@(t, x) F(x) + G(x) * u(x) + U2(x) * u(x) ^ 2 + U3 * u(x) ^ 3, tspan, X0);
u = uCubGarrard; [t3, X3] = ode45(@(t, x) F(x) + G(x) * u(x) + U2(x) * u(x) ^ 2 + U3 * u(x) ^ 3, tspan, X0);

f1 = figure('Position', [64.3333 438 1.6427e+03 420]);
subplot(1, 4, 1); hold on;
plot(t1, X1(:, 1) / (pi / 180))
plot(t2, X2(:, 1) / (pi / 180))
plot(t3, X3(:, 1) / (pi / 180))

subplot(1, 4, 2); hold on;
% plot(t1,X1(:,2)/(pi/180))
% plot(t2,X2(:,2)/(pi/180))
% plot(t3,X3(:,2)/(pi/180))

subplot(1, 4, 3); hold on;
% plot(t1,X1(:,3)/(pi/180))
% plot(t2,X2(:,3)/(pi/180))
% plot(t3,X3(:,3)/(pi/180))

subplot(1, 4, 4); hold on;
% us1 = []; us2 = []; us3 = [];
% for i=1:length(X1),    us1(end+1) = u(X1(i,:).'); end
% for i=1:length(X2),    us2(end+1) = u(X2(i,:).'); end
% for i=1:length(X3),    us3(end+1) = u(X3(i,:).'); end
% plot(t1,us1); plot(t2,us2); plot(t3,us3)

% Second IC
% alpha0 = 30.1;
% alpha0 = 25.628;
alpha0 = 26;
X0 = [pi / 180 * alpha0; 0; 0];

u = uLinGarrard; [t1, X1] = ode45(@(t, x) F(x) + G(x) * u(x) + U2(x) * u(x) ^ 2 + U3 * u(x) ^ 3, tspan, X0);
u = uQuadGarrard; [t2, X2] = ode45(@(t, x) F(x) + G(x) * u(x) + U2(x) * u(x) ^ 2 + U3 * u(x) ^ 3, tspan, X0);
u = uCubGarrard; [t3, X3] = ode45(@(t, x) F(x) + G(x) * u(x) + U2(x) * u(x) ^ 2 + U3 * u(x) ^ 3, tspan, X0);

% figure; hold on;
subplot(1, 4, 1); hold on;
set(gca, 'ColorOrderIndex', 1)

plot(t1, X1(:, 1) / (pi / 180), '--')
plot(t2, X2(:, 1) / (pi / 180), '--')
plot(t3, X3(:, 1) / (pi / 180), '--')

ylim([0 35])
xlim(tspan)
drawnow

subplot(1, 4, 2); hold on;
set(gca, 'ColorOrderIndex', 1)

plot(t1, X1(:, 2) / (pi / 180), '--')
plot(t2, X2(:, 2) / (pi / 180), '--')
plot(t3, X3(:, 2) / (pi / 180), '--')

ylim([-10 5])
xlim(tspan)
drawnow

subplot(1, 4, 3); hold on;
set(gca, 'ColorOrderIndex', 1)

plot(t1, X1(:, 3) / (pi / 180), '--')
plot(t2, X2(:, 3) / (pi / 180), '--')
plot(t3, X3(:, 3) / (pi / 180), '--')

ylim([-15 5])
xlim(tspan)
drawnow

subplot(1, 4, 4); hold on;
us1 = []; us2 = []; us3 = [];
for i = 1:length(X1), us1(end + 1) = u(X1(i, :).'); end
for i = 1:length(X2), us2(end + 1) = u(X2(i, :).'); end
for i = 1:length(X3), us3(end + 1) = u(X3(i, :).'); end
plot(t1, us1 / (pi / 180), '--'); plot(t2, us2 / (pi / 180), '--'); plot(t3, us3 / (pi / 180), '--')

ylim([-10 5])
xlim(tspan)
drawnow

%% Recreate Almubarak Figure 1
eta = 1; % values should be between -\infty and 1.
fprintf('Simulating for eta=%g (gamma=%g)\n', eta, 1 / sqrt(1 - eta))

%  Compute the polynomial approximations to the future energy function
degree = 8;
[w] = pqr(f, g, h, 1 / eta, degree);

tspan = [0, 12];

% First IC
alpha0 = 25;
X0 = [pi / 180 * alpha0; 0; 0];
f3 = figure('Position', [415 47 793.3333 347.3333]);
subplot(1, 2, 1); hold on;
subplot(1, 2, 2); hold on;

legendEntries = {};
for d = 2:2:degree
    u = @(x) (- eta * G(x).' * kronPolyDerivEval(w(1:d), x).' / 2);
    [t, X] = ode45(@(t, x) F(x) + G(x) * u(x) + U2(x) * u(x) ^ 2 + U3 * u(x) ^ 3, tspan, X0);

    subplot(1, 2, 1)
    plot(t, X(:, 1))
    legendEntries{end + 1} = sprintf('order-%i controller ', d - 1);

    subplot(1, 2, 2)
    us = [];
    for i = 1:length(X)
        us(end + 1) = u(X(i, :).');
    end
    plot(t, us)

end

subplot(1, 2, 1)
legend(legendEntries)
ylim([0 0.45])
xlim(tspan)

subplot(1, 2, 2)
legend(legendEntries)
ylim([- .15 0.25])
xlim(tspan)

%% Recreate Krener's Garrard results
tspan = [0, 5];

uLinGarrard = @(x) (- 0.053 * x(1) + 0.5 * x(2) + 0.521 * x(3));
uQuadGarrard = @(x) (- 0.053 * x(1) + 0.5 * x(2) + 0.521 * x(3) + 0.04 * x(1) ^ 2 - 0.048 * x(1) * x(2));
uCubGarrard = @(x) (- 0.053 * x(1) + 0.5 * x(2) + 0.521 * x(3) + 0.04 * x(1) ^ 2 - 0.048 * x(1) * x(2) + 0.374 * x(1) ^ 3 - 0.312 * x(1) ^ 2 * x(2));

%  IC
alpha0 = 22.9;
X0 = [pi / 180 * alpha0; 0; 0];

u = uLinGarrard; [t1, X1] = ode45(@(t, x) F(x) + G(x) * u(x) + U2(x) * u(x) ^ 2 + U3 * u(x) ^ 3, tspan, X0);
u = uQuadGarrard; [t2, X2] = ode45(@(t, x) F(x) + G(x) * u(x) + U2(x) * u(x) ^ 2 + U3 * u(x) ^ 3, tspan, X0);
u = uCubGarrard; [t3, X3] = ode45(@(t, x) F(x) + G(x) * u(x) + U2(x) * u(x) ^ 2 + U3 * u(x) ^ 3, tspan, X0);

f5 = figure('Position', [521.6667 20.3333 736.6667 837.3333]);

% figure; hold on;
subplot(2, 1, 1); hold on;
set(gca, 'ColorOrderIndex', 1)

plot(t1, X1(:, 1) / (pi / 180), '--')
plot(t2, X2(:, 1) / (pi / 180), '--')
plot(t3, X3(:, 1) / (pi / 180), '--')

ylim([0 35])
xlim(tspan)
drawnow

subplot(2, 1, 2); hold on;
us1 = []; us2 = []; us3 = [];
for i = 1:length(X1), us1(end + 1) = u(X1(i, :).'); end
for i = 1:length(X2), us2(end + 1) = u(X2(i, :).'); end
for i = 1:length(X3), us3(end + 1) = u(X3(i, :).'); end
plot(t1, us1 / (pi / 180), '--'); plot(t2, us2 / (pi / 180), '--'); plot(t3, us3 / (pi / 180), '--')

ylim([-10 5])
xlim(tspan)
drawnow

%%
% uLin         = @(x) ( - eta * G(x).' * kronPolyDerivEval(w(1:2), x).'/2);
% uQuad        = @(x) ( - eta * G(x).' * kronPolyDerivEval(w(1:3), x).'/2);
% uCub         = @(x) ( - eta * G(x).' * kronPolyDerivEval(w(1:4), x).'/2);
% set(gca,'ColorOrderIndex',1)
%
% plot(t4,X4(:,1)/(pi/180),':','LineWidth', 2)
% plot(t5,X5(:,1)/(pi/180),':','LineWidth', 2)
% plot(t6,X6(:,1)/(pi/180),':','LineWidth', 2)
% ylim([0 alpha0 + 3])
% xlim(tspan)

end

function [t, y] = ode1(fun, tspan, y0)
% Forward Euler, actually works pretty well
h = 0.01;
t = tspan(1):h:tspan(2);
y = zeros(length(y0), length(t));
y(:, 1) = y0;
for i = 2:length(t)
    y(:, i) = y(:, i - 1) + h * fun(t(i - 1), y(:, i - 1));
end
y = y.';
end

function [t, y] = ode2(odefun, tspan, y0)
% ode2: Implicit Euler ODE solver
h = 0.01;
t = tspan(1):h:tspan(2);
y = zeros(length(y0), length(t));
y(:, 1) = y0;
for i = 2:length(t)
    % Fixed-point iteration for implicit Euler method
    max_iterations = 100;
    y_next = y(:, i - 1);
    for j = 1:max_iterations
        y_next_guess = y(:, i - 1) + h * odefun(t(i), y_next);

        if norm(y_next_guess - y_next) < 1e-12
            break;
        end

        y_next = y_next_guess;
    end
    y(:, i) = y_next;
end
y = y.';
end
