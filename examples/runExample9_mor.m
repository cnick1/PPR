function [] = runExample9_mor(n, r, degree)
%runExample9 Runs the Allen-Cahn example.
%
%   Usage:  [v,w] = runExample9_mor(n,r,degree)
%
%   Inputs: n      - desired state dimension
%           r      - ROM dimension; if r=n, no MOR is performed
%           degree - desired polynomial degree of value function to compute
%   Outputs:
%
%   Background: Based on p34.m from [1].
%
%   Reference: [1] L. N. Trefethen, Spectral methods in MATLAB. Society
%              for Industrial and Applied Mathematics, 2000.
%              doi: 10.1137/1.9780898719598.
%
%   Part of the NLbalancing repository.
%%
if nargin < 3
    degree = 4;
    if nargin < 2
        r = 10;
        if nargin < 1
            n = 33;
        end
    end
end

N = n-1;
fprintf('Running Example 9\n')

%% Construct controller
% Get system
y0 = .5; % Desired interface location

eps = 0.01;
[f, B, ~, D, y] = getSystem9(eps, n-1, y0);
fprintf("Maximum eigenvalue of A is %f; should be zero I think.\n",max(eigs(full(f{1}),n)))

% Reference configuration (@ origin) -> v = v+vref
if isempty(y0)
    vref=zeros(n,1);
else
    vref = tanh((y-y0)/sqrt(2*eps));
end

% Set number of inputs/outputs here

C = B.';

m = 3; p = 3;
B = B(:,linspace(1,n,m+2)); B(:,[1 m+2]) = []; 
C = C(linspace(1,n,p+2),:); C([1 p+2],:) = []; 
Q2 = .1; Q3 = sparse(n^3,1) ; Q4 = sparse(linspace(1,n^4,n),1,4);
q = {[],Q2,Q3,Q4};
R = 1;

fprintf("Computing ppr() solution, n=%i, r=%i, d=%i ... \n",n,r,degree); tic
options.verbose = true; options.r = r; options.eta = 1; options.h = C;
[ValueFun, Gains, options] = ppr(f, B, q, R, degree, options);
fprintf("completed ppr() in %2.2f seconds. \n", toc)

uOpenLoop = @(z) zeros(m,1);
uPPR = @(z) (kronPolyEval(Gains, z));

controllers = {uOpenLoop, uPPR};

for idx = 2%1:2
    u = controllers{idx};
    
    %% Solve PDE by Euler formula and plot results:
    % Construct originial system dynamics
    D2 = D^2; D2([1 n],:) = zeros(2,n); % For boundary conditions
    
    % Initial condition
    v0 = .53*y + .47*sin(-1.5*pi*y);
    % v0 = tanh((y-(-0.125))/sqrt(2*eps*10));
    v = v0;
    
    % Time-stepping
    dt = min([.00001,50*N^(-4)/eps]); t = 0;
    tmax = 20; tplot = 2; nplots = round(tmax/tplot);
    plotgap = round(tplot/dt); dt = tplot/plotgap;
    xx = -1:.025:1; vv = polyval(polyfit(y,v,20),xx);
    plotdata = [vv; zeros(nplots,length(xx))]; tdata = t;
    
    Lagrangian = zeros(length(0:dt:tmax),1);
    
    for i = 1:nplots
        fprintf('%i',i)
        for nn = 1:plotgap
            xbar = v-vref;
            Ux = u(xbar);
            Lagrangian(nn+(i-1)*plotgap) = 1/2*(xbar.'*Q2*xbar + Ux.'*R*Ux + 4*sum(xbar.^4)); % hardcoded v.^4 instead of Q4 for speed
            t = t+dt; v = v + dt*(eps*D2*v + v - v.^3 + B*Ux);    % Euler
        end
        vv = polyval(polyfit(y,v,20),xx);
        plotdata(i+1,:) = vv; tdata = [tdata; t];
    end
    figure, subplot('position',[.1 .4 .8 .5])
    mesh(xx,tdata,plotdata), grid on, axis([-1 1 0 tmax -1.05 1.05]),
    view(-60,55), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
    drawnow
    
    % Compute performance Index (cost)
    performanceIndex = trapz((0:dt:tmax), Lagrangian);
    fprintf("\n\n   The performance index is %f\n\n",performanceIndex)
    
end

% exportgraphics(figure(6),sprintf('plots/example9_mor_n%i_r%i_d%i.pdf',n,r,degree), 'ContentType', 'vector')
% exportgraphics(figure(2),sprintf('plots/example9_mor_n%i_d%i.pdf',n,degree), 'ContentType', 'vector')
% exportgraphics(figure(1),sprintf('plots/example9_mor_openloop_n%i.pdf',n), 'ContentType', 'vector')


end
