% function runExample20()
%runExample20 Runs Vlasov
%
%   Usage:  runExample20()
%
%   Reference: [1]
%
%   Part of the PPR repository.
%%

clear; clc; close all;
epsilon = 0.1;
u = 0; alpha = sqrt(2); Nv = 50; Nx = 20; L = 2*pi; nu = .25; 
NxTotal = 2*Nx+1; v = linspace(-3,3,200); x = linspace(0,L,NxTotal);

[f,g,h] = getSystem20(u,alpha,Nv,Nx,L,nu);
A = f{1}; B = g{1}; K = h{1};
% lyap(A,B*B.')

% FofXU = @(x,phi) (A*x + B*phi);                   % Linearized dynamics (Vlasov equation)
% FofXU = @(x,phi) (A*x + B*phi + g{2}*kron(x,phi));  % Bilinear dynamics   (Vlasov equation)
phi = @(x) sparse((K * x));                                 % Poisson equation

% Define initial condition
x0 = zeros(Nv*NxTotal,1);
x0(Nx) = 0.5 * epsilon / alpha; x0(Nx+2) = 0.5 * epsilon / alpha;

% t = linspace(0,5,5001); % Specify time vector for plotting
t = [0,10]; 

% FofXU(A,B,g{2},x0,phi(x0));

% Simulate using ode solver
% [t, X] = ode15s(@(t, x) FofXU(x,phi(x)),t, x0);
opts_openloop = odeset(OutputFcn=@odeprog);
global T0; T0 = tic;
[t, X] = ode15s(@(t, x) FofXU(A,B,g{2},x,phi(x)),t, x0, opts_openloop);
fprintf("completed in %2.2f seconds. \n", toc(T0))
       
%% Decompose solution into Fourier coefficients 
C = cell(Nv+1,1);

C{1} = ifftshift(X(:,1:NxTotal),2); 
C{1}(:,1) = 1/alpha;
solution = tensorprod(ifft(C{1},[],2),psi(0,alpha,u,v));

for i=1:Nv-1
 C{i+1} = ifftshift(X(:,(1:NxTotal) + i*NxTotal),2); 
 solution = solution + tensorprod(ifft(C{i+1},[],2),psi(i,alpha,u,v));
end
solution = real(solution);

%% Plot a preliminary sanity check
% Electric Field Damping Rate
figure('Position',[65.6667 285 560 420])
title('Electric Field Damping Rate')
semilogy(t,abs(X(:,Nx)),'LineWidth', 2)
hold on 
semilogy(t,abs(X(:,Nx+2)),'--','LineWidth', 2)
semilogy(t,exp(-0.851*t),'k:','LineWidth', 2)



% Solution at initial time
figure('Position',[603 285 560 420])
title('Distribution Function Initial Condition')
mesh(x,v,squeeze(solution(1,:,:)).')
drawnow

% Solution animation
figure('Position',[1117 285 560 420])
for i = 1:1:length(t)
    mesh(x,v,squeeze(solution(i,:,:)).')
    % view(2)
    % zlim([-1e-2 1e-2])
    title(sprintf('Distribution Function vs time; t=%f',t(i)))
    pause(0.01)
end

%% First steps towards a ROM 

fom = ss(A,B,K,[]); 
R = reducespec(fom,'balanced');
% Run algorithm once (optional, recommended for sparse)
R = process(R);
% Select order graphically
view(R)
% Get a reduced model of order 7
rsys = getrom(R,Order=40);



%% Old stuff

% figure
% title('Closed-loop animation')
% for i = 1:10:length(t)
%     % plot(linspace(0,1,220),abs(X(i,:)));
%      % plot(temp(i,:))
%      % ylim([-1e-3,1e-3])
%     % ylim([-0.5 * epsilon / alpha, 0.5 * epsilon / alpha])
% 
%     mesh(x,v,squeeze(solution(i,:,:)).')
%     % view(2)
%     % zlim([-2e-4 2e-4])
%     pause(0.001)
% 
% end
% 
% 
% % Plots just for checking results, not for the paper
% figure
% plotT = t(1:100:end); X = X(1:100:end,:); nplots = length(plotT);
% xx = -1:.025:1; plotdata = [zeros(nplots,length(xx))];
% for i=1:length(plotT)
%     plotdata(i,:) = polyval(polyfit(y,X(i,:),20),xx);
% end
% subplot(2,3,idx)
% mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -1.05 1.05]),
% view(-60,55), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
% title(sprintf("Controller %s  (x0=%2.2f)",controllerNames{idx},x0)); drawnow

% end

function psi_n = psi(n,alpha,u,v)

xi = (v-u)/alpha;
psi_n = (pi*2^n*factorial(n))^(-1/2) * exp(-xi.^2) .* hermiteH(n,xi);
psi_n = psi_n.';
end


function xdot = FofXU(A,B,G2,x,phi)
% FofXU = @(x,phi) (A*x + B*phi + g{2}*kron(x,phi));  % Bilinear dynamics   (Vlasov equation)
% [n,m] = size(B);
% 
% xdot = A*x + B*phi;
% 
% for i=1:n %for each row of G2
%     xdot(i) = xdot(i) + phi.' * reshape(G2(i,:),m,n) * x;
% end

% xdot = (A*x + B*phi + G2*kron(x,phi));
xdot = (A*x + B*phi);

end

function status = odeprog(t, y, flag)
% ODEPROG Custom progress bar for ode solver
% Use with odeset: opts = odeset('OutputFcn',@odeprog);
persistent T1 nSteps lastPct lastUpdateTime
global T0

status = false;

switch flag
    case 'init'
        % Initialize progress bar
        elapsed = toc(T0);
        T1 = t(end);
        nSteps = 50; % Number of blocks in the progress bar
        lastPct = -1;
        lastUpdateTime = 0;
        fprintf(' |%s|  (elapsed: %5i s, remaining: ----- s)', repmat(' ',1,nSteps), round(elapsed));
        
    case ''
        % ODE solver step
        if isempty(t), return; end
        tNow = t(end);
        pct = min(100, max(0, 100 * tNow / T1));
        block = floor(pct / (100/nSteps));
        elapsed = toc(T0);
        eta = (elapsed / max(tNow,eps)) * (T1 - tNow); % avoid divide-by-zero
        needsUpdate = pct - lastPct >= 2 || block == nSteps;
        timeSinceLast = elapsed - lastUpdateTime;
        
        if needsUpdate || timeSinceLast >= 1
            bar = [repmat('-',1,block), repmat(' ',1,nSteps-block)];
            fprintf(repmat('\b',1,93));
            fprintf(' |%s|  (elapsed: %5i s, remaining: %5i s)', bar, round(elapsed), min(round(eta),99999));
            if needsUpdate
                lastPct = pct;
            end
            lastUpdateTime = elapsed;
        end
        
    case 'done'
        % Finalize
        % elapsed = toc(T0);
        % bar = repmat('-',1,nSteps);
        fprintf(repmat('\b',1,93));
        % fprintf(' |%s|  (elapsed: %5i s, remaining:     0 s)\n', bar, round(elapsed));
        clear T0 T1 nSteps lastPct lastUpdateTime
end
end