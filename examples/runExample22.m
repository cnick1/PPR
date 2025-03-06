% function runExample22()
%runExample22 Runs the 4D pendulum on a cart example.
%
%   Usage:  runExample22()
%
%   Inputs:
%
%   References: [1] https://www.mathworks.com/help/control/ref/lti.lqr.html
%               [2] https://www.mathworks.com/help/symbolic/derive-and-simulate-cart-pole-system.html.
%               [3] https://www.mathworks.com/help/mpc/ug/swing-up-control-of-a-pendulum-using-nonlinear-model-predictive-control.html
%
%%
clear; %close all;
global xdot
%% Get model and construct dynamics, initial condition, and time vector
[f, g, h, xdot] = getSystem22(7); % xdot = @(X,u) f{1}*X + g{1}*u;

% Construct original system dynamics
FofXU = @(x,u) xdot(x,u);

% Initial condition
% x0 = [0;0;-pi/3*1.04;0];
% x0 = [0;0;-pi/3*1.06;0];
% x0 = [0;0;-pi/3*1.06;5.12];
% x0 = [0;0;-pi;0]; % ulim = 20;
% x0 = [0;0;-pi;10.975]; ulim = 3;

x0 = [10;0;0;0];

% Time vector
tmax = 10; dt = 0.001; t = 0:dt:tmax-dt; nsteps = length(t);

%% Construct control laws
Q = diag([3 0 3 0]); R = 1; R2 = sparse(1,4^2); R6 = sparse(1,4^6); 
% alpha = pi/10.975; R2(11) = 1; R2(15) = -2*alpha; R2(16) = alpha^2; R2 = 100*R2;
R6(2731) = 0.321;

% PPR/LQR
% [v, GainsPPR] = ppr(f, g, {0,Q,0,0}, {R,0,R2,0,0,0,R6},8);
[v, GainsPPR] = ppr(f, g, {0,Q,0,0}, {R,0,0,0,0,0,0},8);

% MPC controller
nx = 4; nu = 1; nlobj = nlmpc(nx, nx, nu);

nskip = 50;
Ts = dt*nskip; nlobj.Ts = Ts; nlobj.PredictionHorizon = 50; nlobj.ControlHorizon = nlobj.PredictionHorizon;
nlobj.Model.StateFcn = @(x,u,Ts) pendulumDT0(x, u, Ts, nskip);
nlobj.Model.IsContinuousTime = false; nlobj.Model.NumberOfParameters = 1;
nlobj.Weights.OutputVariables = diag(Q).'; nlobj.Weights.ManipulatedVariables = R; nlobj.Weights.ManipulatedVariablesRate = 0.001;
nlobj.OV(1).Min = -30; nlobj.OV(1).Max = 30; nlobj.MV.Min = -300; nlobj.MV.Max = 300;
% validateFcns(nlobj,[0.1;0.2;-pi/2;0.3],0.4,[],{Ts});

nloptions = nlmpcmoveopt; nloptions.Parameters = {Ts};

% Now define control law function handles
uOpenLoop = @(z) zeros(1,1);
uLQR = @(z) (kronPolyEval(GainsPPR, z, 1));
uPPR = @(z) (kronPolyEval(GainsPPR, z));
uMPC = @(x,uprev,nloptions) nlmpcmove(nlobj,x,uprev,[0 0 0 0],[],nloptions); 

controllers = {uOpenLoop, uLQR, uPPR, uMPC};
controllerNames = {'Uncontrolled', 'LQR', 'PPR', 'MPC'};

%% Simulate closed-loop systems
fig1 = figure('Position',[106.3333 308.3333 862.0000 544.0001]); fig2 = figure('Position',[228.3333 30.3333 674 261.3333]);
hbar = waitbar(0,'Simulation Progress');
for idx = 1:3
    u = controllers{idx}; 
    x=x0; X = []; Ux = []; Lagrangian = zeros(size(t)); 
    fprintf("Performing closed-loop simulation with %s controller ... ",controllerNames{idx}); tic
    for ct=1:nsteps
        % Compute optimal control moves
        uprev = u(x); 
        if exist('ulim','var'); if uprev > ulim; uprev = ulim; elseif uprev < -ulim; uprev=-ulim; end; end % actuator limits sometimes help PPR controller
        X = [X; x.']; Ux = [Ux; uprev];
        Lagrangian(ct) = 1/2*(x.'*Q*x + uprev.'*R*uprev); 

        % Time step
        x = x + dt*xdot(x,uprev);  

        waitbar(ct/nsteps,hbar);
    end
    fprintf("completed in %2.2f seconds. \n", toc)
    performanceIndex(idx) = trapz(t, Lagrangian); drawnow
    figure(fig1); for i=1:4; subplot(2,2,i); hold on; plot(t,X(:,i)); end
    figure(fig2); hold on; plot(t,Ux);
end
%% Add plot labels, limits, etc.
figure(fig1); subplot(2,2,1); xlabel('time'); ylabel('z'); title('cart position'); ylim([-4 4])
subplot(2,2,2); xlabel('time'); ylabel('zdot'); title('cart velocity'); ylim([-5 7])
subplot(2,2,3); xlabel('time'); ylabel('theta'); title('pendulum angle'); ylim([-2*pi 2*pi])
subplot(2,2,4); xlabel('time'); ylabel('thetadot'); title('pendulum velocity'); ylim([-10 10])
figure(fig2); xlabel('time'); ylabel('u'); title('control signal'); ylim([-3 3])
drawnow


%%
x=x0.'; X = []; Ux = []; Lagrangian = []; uprev = 0; 
fprintf("Performing closed-loop simulation with MPC controller with nskip=%i, time horizon=%3.2f ... ",nskip,nlobj.PredictionHorizon*Ts); tic
for ct=1:nsteps/nskip
    % Compute optimal control moves
    [uprev,nloptions] = uMPC(x,uprev,nloptions);
    
    % Time step (using one MPC control for multiple steps due to computational burden)
    for ct2=1:nskip
        X = [X; x]; Ux = [Ux; uprev];
        Lagrangian = [Lagrangian; 1/2*(x*Q*x.' + uprev.'*R*uprev)]; 
        x = x + dt*xdot(x.',uprev).';
    end

    waitbar(ct/nsteps*nskip,hbar);
end
fprintf("completed in %2.2f seconds. \n", toc)
performanceIndex(4) = trapz(t, Lagrangian); drawnow
figure(fig1); for i=1:4; subplot(2,2,i); hold on; plot(t,X(:,i)); end
figure(fig2); hold on; plot(t,Ux);

%% Add plot labels, limits, etc.
figure(fig1); subplot(2,2,1); xlabel('time'); ylabel('z'); title('cart position'); ylim([-1 3])
subplot(2,2,2); xlabel('time'); ylabel('zdot'); title('cart velocity'); ylim([-6.5 6.5])
subplot(2,2,3); xlabel('time'); ylabel('theta'); title('pendulum angle'); ylim([-2*pi 2*pi])
subplot(2,2,4); xlabel('time'); ylabel('thetadot'); title('pendulum velocity'); ylim([-10 10])
legend(controllerNames)

figure(fig2); xlabel('time'); ylabel('u'); title('control signal u(t)'); ylim([-30 40]); xlim([0 4]);
drawnow
legend(controllerNames)
close(hbar)

%% Print cost table
fprintf('\n# Table I Data (pendulum on a cart)\n');
fprintf('# Control costs for different controllers\n');
fprintf("      %s    &    %s    &    %s    &    %s  \n",controllerNames{1},controllerNames{2},controllerNames{3},controllerNames{4})
fprintf("       %9.3f      &  %6.3f   &  %6.3f   &  %6.3f       \n",performanceIndex)

function xk1 = pendulumDT0(xk, uk, Ts, nskip)
%% Discrete-time nonlinear dynamic model of a pendulum on a cart at time k
global xdot
% Repeat application of Euler method sampled at Ts/nskip.
delta = Ts/nskip;
xk1 = xk;
for ct=1:nskip
    xk1 = xk1 + delta*xdot(xk1,uk);
end
end

% end
