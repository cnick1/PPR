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
u = 0; alpha = sqrt(2); Nv = 50; Nx = 20; L = 2*pi; nu = .5; 
NxTotal = 2*Nx+1; v = linspace(-3,3,200); x = linspace(0,L,NxTotal);

[f,g,h] = getSystem20(u,alpha,Nv,Nx,L,nu);
A = f{1}; B = g{1}; K = h{1};
% lyap(A,B*B.')


FofXU = @(x,phi) (A*x + B*phi);                   % Linearized dynamics (Vlasov equation)
% FofXU = @(x,phi) (A*x + B*phi + g{2}*kron(x,phi));  % Bilinear dynamics   (Vlasov equation)
phi = @(x) (K * x);                                 % Poisson equation

% Define initial condition
x0 = zeros(Nv*NxTotal,1);
x0(Nx) = 0.5 * epsilon / alpha; x0(Nx+2) = 0.5 * epsilon / alpha;

% t = linspace(0,5,5001); % Specify time vector for plotting
t = [0,30]; 

% Simulate using ode solver
[t, X] = ode15s(@(t, x) FofXU(x,phi(x)),t, x0);

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
title('Distribution Function vs time')
for i = 1:10:length(t)
    mesh(x,v,squeeze(solution(i,:,:)).')
    % view(2)
    % zlim([-1e-2 1e-2])
    pause(0.001)
end

%% Old stuff

figure
title('Closed-loop animation')
for i = 1:10:length(t)
    % plot(linspace(0,1,220),abs(X(i,:)));
     % plot(temp(i,:))
     % ylim([-1e-3,1e-3])
    % ylim([-0.5 * epsilon / alpha, 0.5 * epsilon / alpha])

    mesh(x,v,squeeze(solution(i,:,:)).')
    % view(2)
    % zlim([-2e-4 2e-4])
    pause(0.001)

end


% Plots just for checking results, not for the paper
figure
plotT = t(1:100:end); X = X(1:100:end,:); nplots = length(plotT);
xx = -1:.025:1; plotdata = [zeros(nplots,length(xx))];
for i=1:length(plotT)
    plotdata(i,:) = polyval(polyfit(y,X(i,:),20),xx);
end
subplot(2,3,idx)
mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -1.05 1.05]),
view(-60,55), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
title(sprintf("Controller %s  (x0=%2.2f)",controllerNames{idx},x0)); drawnow

% end

function psi_n = psi(n,alpha,u,v)

xi = (v-u)/alpha;
psi_n = (pi*2^n*factorial(n))^(-1/2) * exp(-xi.^2) .* hermiteH(n,xi);
psi_n = psi_n.';
end
