function runExample27()
%runExample27 Runs the 1D unstable heat equation FEM control example.
%
%   Usage:  runExample27()
%
%   Description: The model, inspired by [1], describes an unstable heat
%   equation problem. The model takes the form
%
%     uₜ(x,t)  = ε uₓₓ(x,t) + λ u(x,t) + μ u(x,t)³
%     uₓ(0,t) = 0                     (Neumann BC)
%     uₓ(1,t) = u(t)                  (Neumann boundary control input)
%
%   where one end is insulated and one end is subject to Neumann boundary
%   control. The FEM model can be written (after multiplying by M⁻¹) as
%
%       ẋ = A x + F₃(x⊗x⊗x) + B u
%       y = C x
%
%   for which we can compute a controller u(x) = K(x) using PPR.
%
%   Reference: [1] D. M. Boskovic, M. Krstic, and W. Liu, "Boundary control
%              of an unstable heat equation via measurement of
%              domain-averaged temperature,” IEEE Transactions on Automatic
%              Control, vol. 46, no. 12, pp. 2022–2028, 2001
%              [2] S. A. Ragab and H. E. Fayed, Introduction to finite
%              element analysis for engineers. Taylor & Francis Group, 2017
%
%   Part of the PPR repository.
%%
fprintf('Running Example 27\n')

% Get dynamics
n = 40;
[f, g, ~, xg] = getSystem27(n-1,.75,1,-1);
F = @(x) kronPolyEval(f, x); G = @(x) g{1};
% F = @(x) f{1}*x; G = @(x) g{1};

% Get value function/controller
q = 0.5; R = 1; degree = 4;
fprintf(" Computing ppr() solution, n=%i, d=%i ... ",n,4); tic
[~, K] = ppr(f, g, q, R, degree);
fprintf(" completed in %2.2f seconds. \n", toc)

uLQR = @(x) kronPolyEval(K, x, 1);
uPPR = @(x) kronPolyEval(K, x, degree-1);

% To compare with back-stepping controller from [1]
% a = 0.95*pi/2; alpha = 3;
% uBS = @(x) -(alpha + a*tan(a))*x(end) - (alpha*a*tan(a) + a^2/cos(a)^2)*trapz(xg,x);

%% Simulate closed-loop system
% alpha = 4; d = -8;
% alpha = 12; d = -1;
% b = [1 1;2 3]\[-1-d;-alpha*d]; c=b(2);b=b(1);
% x0 = 1 + b*xg.^2 + c*xg.^3 + d*xg.^alpha;
% figure;plot(xg,x0)
% [b,c,d]
% x0 = 1 +b*xg.^2 + c* xg.^3 + d*xg.^alpha;


% x0 = 1 - 11*xg.^2 + 18* xg.^3 - 8*xg.^4;
% x0 = 1 - 15*xg.^2 + 26* xg.^3 - 12*xg.^4;
% x0 = 1 - 13.5*xg.^2 + 14* xg.^3 - 1.5*xg.^10;
x0 = 1 - 12*xg.^2 + 12* xg.^3 - 1*xg.^12;
% x0 = 1.1 + 0.*xg;
tmax = 2; t = 0:0.02:tmax; % specify for plotting

[~, XUNC] = ode15s(@(t, x) F(x)                 , t, x0);
[~, XLQR] = ode15s(@(t, x) F(x) + G(x) * uLQR(x), t, x0);
[t, XPPR] = ode15s(@(t, x) F(x) + G(x) * uPPR(x), t, x0);


%% Plot solution
figure('Position',[501 19 632 838])
% figure('Position',[474.3333 340.3333 925.3333 300.6667]);
subplot(3,1,1)
mesh(xg,t(1:size(XUNC,1)),XUNC);
grid on, axis([0 1 0 tmax -1.125 1.125]), view(145,35), colormap([0 0 0]);
xlabel x, ylabel t, zlabel u(x,t), title Uncontrolled; drawnow

% figure('Position',[474.3333 340.3333 925.3333 300.6667]);
subplot(3,1,2)
mesh(xg,t(1:size(XLQR,1)),XLQR);
grid on, axis([0 1 0 tmax -1.125 1.125]), view(145,35), colormap([0 0 0]);
xlabel x, ylabel t, zlabel u(x,t), title LQR; drawnow


subplot(3,1,3)
mesh(xg,t,XPPR);
grid on, axis([0 1 0 tmax -1.125 1.125]), view(145,35), colormap([0 0 0]);
xlabel x, ylabel t, zlabel u(x,t), title PPR; drawnow


% figure
% % Animate
% for i=1:10:length(t)
%     plot(xg,X(i,:));
%     xlim([0 1]); ylim([-1 1]); drawnow
% end



end



