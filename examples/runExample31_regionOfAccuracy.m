function runExample31_regionOfAccuracy()
%runExample31_regionOfAccuracy Runs Fokker-Planck example to compare optimal
% control computation. This function plots error vs region size comparing
% degree of approximation for a polynomial function.
%
%   Usage:  runExample31_regionOfAccuracy()
%
%   Part of the PPR repository.

[f, g, h, x0] = getSystem31(false, '2');
q = h{1}.'*h{1}; R = 1e-4;
n = length(x0);
FofXU = @(x,u) kronPolyEval(f,x) + kronPolyEval(g,x,scenario='G(x)') * u;

%% Compute true PPR Controllers (only up to degree 3)
% Get value function/controller
fprintf(" Computing ppr() solution, n=%i d=%i ... ",n,4); tic
[v, K] = ppr(f, g, q, R, 4);
fprintf(" completed in %2.2f seconds. \n", toc)
uLQR = @(x) kronPolyEval(K, x, degree=1);
uPPR2 = @(x) kronPolyEval(K, x, degree=2);
uPPR3 = @(x) kronPolyEval(K, x, degree=3);

% Set up time vector for simulations
tmax = 20; t = linspace(0,tmax,2001);
options_ode = odeset('RelTol',1e-8);

% Set up performance index (cost)
runningCost = @(x,u) 1/2 * sum((x*q) .* x, 2) + 1/2 * R * u.^2;

% Set up ||x|| and control cost vectors
xNorms = 4:0.1:5;
OptimalCost = zeros(size(xNorms));
PPR1Cost = zeros(size(xNorms));
PPR2Cost = zeros(size(xNorms));
PPR3Cost = zeros(size(xNorms));

uOL = zeros(length(t),1); % Use uOL as initial guess for optimal control; initialize to zero

for i = 1:length(xNorms)
    
    %% Optimal control computed as open-loop signal using gradient descent via a
    % two-point boundary value/Pontryagin-like first-order optimality approach
    fprintf(" - Simulating optimal closed-loop dynamics via gradient descent, x0=%.1f ... ",xNorms(i)); T0 = tic;
    sigma= 0.05; theta= 0.7; tolerance= 1e-9;
    [xOL, uOL] = optimize(uOL.',f{1},g{2},g{1},h{1},R,xNorms(i)*x0,t, ...
        odeset('RelTol',1e-13,'MaxOrder',2),sigma,theta,tolerance);
    fprintf(" completed in %2.2f seconds. \n", toc(T0))
    
    % xOL = zeros(length(t),n); uOL = 0*t.';
    
    %% Simulate PPR closed-loop system
    fprintf(" - Simulating LQR closed-loop dynamics ... "); T0 = tic;
    [~, XLQR] = ode15s(@(t, x) FofXU(x, uLQR(x)), t, xNorms(i)*x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))
    fprintf(" - Simulating PPR2 closed-loop dynamics ... "); T0 = tic;
    [~, XPPR2] = ode15s(@(t, x) FofXU(x, uPPR2(x)), t, xNorms(i)*x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))
    fprintf(" - Simulating PPR3 closed-loop dynamics ... "); T0 = tic;
    [~, XPPR3] = ode15s(@(t, x) FofXU(x, uPPR3(x)), t, xNorms(i)*x0, options_ode); fprintf("completed in %2.2f seconds. \n", toc(T0))
    
    %% Compute control signals and integrated costs
    UxLQR = zeros(length(t),1);
    UxPPR2 = zeros(length(t),1);
    UxPPR3 = zeros(length(t),1);
    for j=1:length(t)
        UxLQR(j,1) = uLQR(XLQR(j,:).');
        UxPPR2(j,1) = uPPR2(XPPR2(j,:).');
        UxPPR3(j,1) = uPPR3(XPPR3(j,:).');
    end
    
    OptimalCost(i) = trapz(t, runningCost(xOL,uOL));
    PPR1Cost(i) = trapz(t, runningCost(XLQR,UxLQR));
    PPR2Cost(i) = trapz(t, runningCost(XPPR2,UxPPR2));
    PPR3Cost(i) = trapz(t, runningCost(XPPR3,UxPPR3));
end

%% Plot J(x0,u) - J(x0,u*) vs. ||x0|| for each controller
figure; hold on;
plot(xNorms*norm(x0),PPR1Cost-OptimalCost,'DisplayName','$u_1(t)$')
plot(xNorms*norm(x0),PPR2Cost-OptimalCost,'DisplayName','$u_2(t)$')
plot(xNorms*norm(x0),PPR3Cost-OptimalCost,'DisplayName','$u_3(t)$')

%% Export data
fileName = sprintf('plots/example31_regionOfAccuracy.dat');
fprintf("Writing data to " + fileName + '\n')
fileID = fopen(fileName, 'w');

%print the header
fprintf(fileID, ['Initial Condition Norm      & PPR1 Cost Gap    & PPR2 Cost Gap    & PPR3 Cost Gap    ' ...
    '& NormPPR1 Cost Gap & NormPPR2 Cost Gap    & NormPPR3 Cost Gap    ' ...
    '& Optimal Cost   & PPR1 Cost        & PPR2 Cost        & PPR3 Cost        \n']);
for ii = 1:length(xNorms)
    fprintf(fileID, '%26.6e      & ', xNorms(ii)*norm(x0));
    fprintf(fileID, '%8.6e     & %8.6e     & %8.6e     & ', ...
        PPR1Cost(ii)-OptimalCost(ii), PPR2Cost(ii)-OptimalCost(ii), PPR3Cost(ii)-OptimalCost(ii));
    fprintf(fileID, '%8.6e     & %8.6e     & %8.6e     & ', ...
        (PPR1Cost(ii)-OptimalCost(ii))/OptimalCost(ii), (PPR2Cost(ii)-OptimalCost(ii))/OptimalCost(ii), (PPR3Cost(ii)-OptimalCost(ii))/OptimalCost(ii));
    fprintf(fileID, '%8.6e     & %8.6e     & %8.6e     & %8.6e     \n', ...
        OptimalCost(ii), PPR1Cost(ii), PPR2Cost(ii), PPR3Cost(ii));
end
fclose(fileID);
type(fileName);

end

%% The code below was provided by Breiten et al.
function [x, u] = optimize(uInit,A,N,B,C,gamma,xInit,tsim,options,sigma,beta,tolerance)
figure;
u= uInit;
tsimBack= fliplr(tsim);
[~,x] = ode15s(@state_equation,tsim,xInit,options,A,N,B,u,tsim);
x= x';
[~,p] = ode15s(@costate_equation,tsimBack,zeros(length(xInit),1),options,A,N,C,x,u,tsim);
p= p';
p = fliplr(p);
gradU = gamma*u + B'*p + ones(1,size(x,1))*(p.*(N*x));
error= L2norm_gradient(gradU, tsim);
test= (L2norm_gradient(gradU, tsim)< tolerance);

costSaved= costFunction(u,A,N,B,C,gamma,xInit,tsim,options);
k= 0;

while (~test)&&(k<85)
    k= k+1;
    step_length = compute_step_length(u,gradU,A,N,B,C,gamma,xInit,tsim,options,sigma,beta);
    u= u - step_length*gradU;
    [~,x] = ode15s(@state_equation,tsim,xInit,options,A,N,B,u,tsim);
    x= x';
    [~,p] = ode15s(@costate_equation,tsimBack,zeros(length(xInit),1),options,A,N,C,x,u,tsim);
    p= p';
    p = fliplr(p);
    gradU = gamma*u + B'*p + ones(1,size(x,1))*(p.*(N*x));
    error= L2norm_gradient(gradU, tsim); semilogy(k,error,'*'); hold on; grid on;
    ylabel('Error'); xlabel('Iteration'); title('Gradient Descent Convergence')
    drawnow
    test= (L2norm_gradient(gradU, tsim)< tolerance);
end
close
costNew= costFunction(u,A,N,B,C,gamma,xInit,tsim,options);
costSaved-costNew;

x=x.';
u=u.';

end

function step_length = compute_step_length(u,gradU,A,N,B,C,gamma,xInit,tsim,options,sigma,beta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

currentCost = costFunction(u,A,N,B,C,gamma,xInit,tsim,options);
slope = L2norm_gradient(gradU, tsim);
step_length= 500;
newCost= costFunction(u-step_length*gradU,A,N,B,C,gamma,xInit,tsim,options);
test= (newCost <= currentCost - sigma*step_length*slope);

while (~test)
    step_length= beta*step_length;
    newCost= costFunction(u-step_length*gradU,A,N,B,C,gamma,xInit,tsim,options);
    test= (newCost <= currentCost - sigma*step_length*slope);
end

variation= newCost-currentCost;


end

function dx = state_equation(t,x,A,N,B,uvalues,tsim)

u = interp1(tsim,uvalues',t);
dx = A*x + N*x*u +B*u;
%t
end

function dp = costate_equation(t,p,A,N,C,xvalues,uvalues,tsim)

x = interp1(tsim,xvalues',t)';
u = interp1(tsim,uvalues',t);
dp = -((C'*C)*x+ (A'+(u*N'))*p);
%t
end

function slope = L2norm_gradient(gradU, tsim)

fun= @(t) (interp1(tsim,gradU',t)).^2;
slope = integral(fun,0,tsim(end));

end

function cost = costFunction(u,A,N,B,C,gamma,xInit,tsim,options)

[~,x] = ode23(@state_equation,tsim,xInit,options,A,N,B,u,tsim);
x= x';
fun= @(t) integral_cost(t,C,gamma,x,u,tsim);
cost= integral(fun,0, tsim(end));

end

function [ TA ] = mat_mult( T,d,n,mu,A )
%MAT_MULT Summary of this function goes here
%   Detailed explanation goes here

m = size(A,1);
TA = reshape(T,n*ones(1,d));
TA = A*reshape(permute(TA,[mu,setdiff(1:d,mu)]),n,n^(d-1));
TA = reshape(TA,[m,n*ones(1,d-1)]);
TA = ipermute(TA,[mu,setdiff(1:d,mu)]);

end


function [ TA ] = mat_mult_all( T,A,d )
%MAT_MULT Summary of this function goes here
%   Detailed explanation goes here

m = size(A,1);
n = size(A,2);

TA = T;
for i = 1:d
    TA = A*reshape(TA,n,numel(TA)/n);
    TA = reshape(TA,[m*ones(1,i),n*ones(1,d-i)]);
    TA = permute(TA,[2:d,1]);
end

end

function z = integral_cost(t,C,gamma,xvalues,uvalues,tsim)

u = interp1(tsim,uvalues',t);
x = interp1(tsim,xvalues',t)';

z= 0.5*(ones(1,size(C,1))*((C*x).^2) + gamma*u.^2);

end

