function [] = runExample9(n, degree)
%runExample9 Runs the Allen-Cahn example.
%
%   Usage:  [v,w] = runExample9()
%
%   Inputs: n      - desired state dimension
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
if nargin < 2 
    degree = 4;
    if nargin < 1
        n = 32;
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

B = B(:,linspace(1,n,5)); B(:,[1 5]) = []; m = size(B,2);
Q2 = .1; Q3 = sparse(n^3,1) ; Q4 = sparse(linspace(1,n^4,n),1,4);
q = {[],Q2,Q3,Q4};
R = 1;

fprintf("Computing ppr() solution, n=%i, d=%i ... \n",n,degree)
[ValueFun, Gains] = ppr(f, B, q, R, degree, true);
fprintf("completed.\n")

uOpenLoop = @(z) zeros(m,1);
uPPR = @(z) (kronPolyEval(Gains, z));

controllers = {uOpenLoop, uPPR};

for idx = 1:2
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
    tmax = 100; tplot = 2; nplots = round(tmax/tplot);
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
    mesh(xx,tdata,plotdata), grid on, axis([-1 1 0 tmax -1 1]),
    view(-60,55), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
    drawnow

    % Compute performance Index (cost)
    performanceIndex = trapz((0:dt:tmax), Lagrangian);
    fprintf("\n\n   The performance index is %f\n\n",performanceIndex)

end

exportgraphics(gcf,sprintf('plots/example9_n%i_d%i.pdf',n,degree), 'ContentType', 'vector')
close
exportgraphics(gcf,sprintf('plots/example9_openloop_n%i.pdf',n), 'ContentType', 'vector')
close


end
