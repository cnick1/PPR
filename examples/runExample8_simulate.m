function [w] = runExample8_simulate(exportData, x0, varargin)
%runExample8_simulate Runs the finite element heat equation example to demonstrate
%            convergence and scalability.
%
%   Usage:  [v,w] = runExample8_simulate(plotEnergy,plotBalancing,balancingDegree,
%                               numGTermsModel, numGTermsApprox, exportData)
%
%   Inputs:
%       exportData      - Boolean variable to determine if
%                         plots/data are exported
%
%   Outputs:
%       v,w             - Coefficients of the past and future energy
%                         function approximations, respectively
%
%   The value of eta is set below.
%
%   Part of the NLbalancing repository.
%%

if nargin < 2
    if nargin < 1
        exportData = false;
    end
    x0 = 1e-5;
end

%% Open Loop Simulations
fprintf('Simulating Example 8, open loop\n')

numEls = [16, 32];
for numEl = numEls
    %% Get dynamics
    [~,~,~,~, f, g, h] = getSystem8(numEl);
    %     f = {f{1}};
    
    % Construct initial condition (from Mark Embree's talk)
    L = 30; x = linspace(0,L,numEl+1).';
    initialCondition = x0 * x .* (x-L).*(x-L/2);
        
    %% Simulate model with ode45
    % Remove end nodes since they are fixed
    X0 = initialCondition(2:end-1);
    tspan  = [0, 120];
    
    tic
    [t1, X1] = ode45(@(t,x) kronPolyEval(f, x), tspan, X0);
    toc
     
    n = length(f{1});
    ff = {}; ff{1} = tensor(f{1}); ff{2} = tensor(reshape(full(f{2}),n,n,[])); ff{3} = tensor(reshape(full(f{3}),n,n,n,[]));
    tic
    [t2, X2] = ode45(@(t,x) tensPolyEval(ff, x), tspan, X0);
    toc
    
    %% Plot figures
    % Add first and last node zero BCs back
    X1 = [zeros(length(t1),1), X1, zeros(length(t1),1)];
    
    figure 
    [T,X] = meshgrid(t1,x);
    surf(T,X,X1.','FaceColor','interp','EdgeColor','none')
    set(gca,'Ydir','reverse')
    
    %%Plot the mesh lines
    xnumlines = 20; xspacing = round(length(x)/xnumlines); % 10 lines
    ynumlines = 20; yspacing = round(length(t1)/ynumlines); % 10 partitions
    
    hold on
    for i = 1:yspacing:length(t1); plot3(t1(i)*ones(size(x)),x,X1(i,:),'-k'); end
    for i = 1:xspacing:length(x); plot3(t1,x(i)*ones(size(t1)),X1(:,i),'-k'); end
    
    drawnow
    
    % Draw a similar figure using waterfall command
    figure
    waterfall(X(:,1:20:end).',T(:,1:20:end).',X1(1:20:end,:)), view(10,70), colormap(1e-6*[1 1 1]);
    drawnow
    
    
end

%% Closed loop simulations
fprintf('Simulating Example 8, closed loop\n')

numEls = 32;
[~,~,~,~, f, g, h] = getSystem8(numEl);
    
degree = 4;
eta = 1;

[w] = pqr(f, g, h2q(h), eta, degree, true);


end
