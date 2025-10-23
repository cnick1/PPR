%  This script runs the Burgers equation control example.
%
%   Reference: [1] B. Kramer, S. Gugercin, J. Borggaard, and L. Balicki, “Nonlinear
%               balanced truncation: Part 1—computing energy functions,” arXiv,
%               Dec. 2022. doi: 10.48550/ARXIV.2209.07645
%              [2] B. Kramer, S. Gugercin, and J. Borggaard, “Nonlinear balanced
%               truncation: Part 2—model reduction on manifolds,” arXiv, Feb. 2023.
%               doi: 10.48550/ARXIV.2302.02036
%              [3] J. Borggaard and L. Zietsman, “On approximating polynomial-
%               -quadratic regulator problems,” IFAC-PapersOnLine, vol. 54, no. 9,
%               pp. 329–334, 2021, doi: 10.1016/j.ifacol.2021.06.090
%
%  Part of the PPR repositories.
%%
clear; clc; 
% close all; 

n = 52; m = 4; 

epsilon = 0.0002;
alpha = 0.02; % alpha = 1.0, not controllable at n=8

[f, g, ~, zInit,M] = getSystem3(n, m, 1, epsilon, alpha);
% plot(g{1})

y = linspace(0,1,n+2); y = y(2:end-1);
% x0 = .01*exp(-(y-0.5).^2/2/0.005); % Gaussian pulse centered at 0.5
% x0 = .0075*exp(-(y-.6).^2/2/0.0025) - .0075*exp(-(y-.4).^2/2/0.0025); % Gaussian pulse centered at 0.5
% x0 = .01*exp(-(y-.55).^2/2/0.001) - .01*exp(-(y-.45).^2/2/0.001); % Gaussian pulse centered at 0.5
% x0 = x0.*1e-2;
% x0 = sqrt(M)*x0.';

% zInit = 0.5*sin(2*pi*x)^2 on (0,0.5) and 0 otherwise
x0 = 0.1*sin(2*pi*y).^2./(2*2*pi);
x0(n/2:end) = 0;

% x0 = 0.01*zInit;

% Q = full(M); 
Q = 1; R = 1;

%% Construct controllers
% Open-loop (uncontrolled) controller
uUnc = @(z) zeros(m,1);

% PPR Controller
fprintf("Computing ppr() solution, n=%i, d=%i ... ",n,4); tic
options = struct; options.verbose = true; 
[~, GainsPPR1] = ppr(f, g, Q, R, 4, options); fprintf("completed in %2.2f seconds. \n", toc)
uLQR = @(z) (kronPolyEval(GainsPPR1, z, 1));
uPPR1 = @(z) (kronPolyEval(GainsPPR1, z));

% [~, GainsPPR2] = ppr(f, g, {0,Q,0,0}, {R,0}, 4, options);
uPPR2 = @(z) (kronPolyEval(GainsPPR1, z));

fprintf("Computing ppr() solution, n=%i, d=%i ... ",n,4); tic
options = struct; options.verbose = true; options.reducedDimension = 14;
[~, GainsPPR3] = ppr(f, g, Q, R, 4, options); fprintf("completed in %2.2f seconds. \n", toc)
uPPR3 = @(z) (kronPolyEval(GainsPPR3, z));

%% Simulate closed-loop systems
% Construct original system dynamics
FofXU = @(v,u) (f{1}*v + f{2}*kron(v,v) + g{1}*u);

t = [0,50]; % opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-20);

j=1;
% Simulate using ode solver (more efficient than forward euler) and
% compute performance indexes
[tUnc, Xunc] = ode15s(@(t, v) FofXU(v,uUnc(v)),t, x0);
costUnc(j) = trapz(tUnc, sum((Xunc.^2).*diag(Q).', 2));

[tLQR, XLQR] = ode15s(@(t, v) FofXU(v,uLQR(v)),t, x0);
for i=1:length(tLQR); UxLQR(i,:) = uLQR(XLQR(i,:).'); end
costLQR(j) = trapz(tLQR, sum((XLQR.^2).*diag(Q).', 2) + sum(R*UxLQR.^2,2));

[tPPR1, XPPR1] = ode15s(@(t, v) FofXU(v,uPPR1(v)),t, x0);
for i=1:length(tPPR1); UxPPR1(i,:) = uPPR1(XPPR1(i,:).'); end
costPPR1(j) = trapz(tPPR1, sum((XPPR1.^2).*diag(Q).', 2) + sum(R*UxPPR1.^2,2));

[tPPR2, XPPR2] = ode15s(@(t, v) FofXU(v,uPPR2(v)),t, x0);
for i=1:length(tPPR2); UxPPR2(i,:) = uPPR2(XPPR2(i,:).'); end
costPPR2(j) = trapz(tPPR2, sum((XPPR2.^2).*diag(Q).', 2) + sum(R*UxPPR2.^2,2));

[tPPR3, XPPR3] = ode15s(@(t, v) FofXU(v,uPPR3(v)),t, x0);
for i=1:length(tPPR3); UxPPR3(i,:) = uPPR3(XPPR3(i,:).'); end
costPPR3(j) = trapz(tPPR3, sum((XPPR3.^2).*diag(Q).', 2) + sum(R*UxPPR3.^2,2));

% switch plots
%     case 'animation'
% fig2 = figure('Position',[664.3333 1.5323e+03 1.0133e+03 557.3333]);
% plot(tUnc,tUnc*0)
% hold on;
% plot(tLQR,UxLQR)
% plot(tPPR1,UxPPR1)
% plot(tPPR2,UxPPR2)
% plot(tPPR3,UxPPR3)
% 
% legend('Uncontrolled','LQR','PPR1','PPR2','PPR3')
% xlim([0 5]);ylim([-5 1])
% drawnow

% Xunc = Xunc/sqrt(M);

figure
title('Closed-loop animation')
for j2 = 1:length(tUnc)
    tt = tUnc(j2);
    hold off
    [~,i] = min(abs(tUnc-tt));
    plot(y,Xunc(i,:));
    hold on
    [~,i] = min(abs(tLQR-tt));
    plot(y,XLQR(i,:));
    [~,i] = min(abs(tPPR1-tt));
    plot(y,XPPR1(i,:));
    [~,i] = min(abs(tPPR2-tt));
    plot(y,XPPR2(i,:));
    [~,i] = min(abs(tPPR3-tt));
    plot(y,XPPR3(i,:));
    ylim([-0.01 0.01]); xlim([0 1])
    legend('Uncontrolled','LQR','PPR1','PPR2','PPR3')
    pause(0.01)
end
%     case 'figs'
%         fig2 = figure;
%         plot(tUnc,tUnc*0)
%         hold on;
%         plot(tPPR2,UxPPR2)
%         plot(tLQR,UxLQR)
%         plot(tPPR1,UxPPR1)
%         % plot(tPPR3,UxPPR3)
%         % plot(tPPR4,UxPPR4)
%         legend('Uncontrolled','PPR2','PPR','PPR reduced')%,'PPR3','TT-HJB')
%
%         figure('Position',[113 277.6667 1.5807e+03 482.6667])
%         % Plots just for checking results, not for the paper
%         plotIndices = round(linspace(0,1,51).^3*5000+1);
%         % plotIndices = linspace(1,5001,51);
%         plotT = t(1:100:end); nplots = length(plotT);
%         xx = -1:.025:1; plotdata = [zeros(nplots,length(xx))];
%         X = Xunc(plotIndices,:);
%         for i=1:length(plotT)
%             plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
%         end
%         subplot(2,3,1)
%         mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -2.05 3.05]),
%         view(-60,55), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
%         title("Open-loop"); drawnow
%
%         X = XPPR2(plotIndices,:);
%         for i=1:length(plotT)
%             plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
%         end
%         subplot(2,3,2)
%         mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -2.05 3.05]),
%         view(-60,55), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
%         title("PPR2"); drawnow
%
%         X = XLQR(plotIndices,:);
%         for i=1:length(plotT)
%             plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
%         end
%         subplot(2,3,3)
%         mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -2.05 3.05]),
%         view(-60,55), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
%         title("PPR"); drawnow
%
%         X = XPPR1(plotIndices,:);
%         for i=1:length(plotT)
%             plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
%         end
%         subplot(2,3,4)
%         mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -2.05 3.05]),
%         view(-60,55), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
%         title("PPR reduced"); drawnow
%
%         X = XPPR3(plotIndices,:);
%         for i=1:length(plotT)
%             plotdata(i,:) = polyval(polyfit(y,X(i,:),8),xx);
%         end
%         subplot(2,3,5)
%         mesh(xx,plotT,plotdata), grid on, axis([-1 1 0 tmax -2.05 3.05]),
%         view(-60,55), colormap([0 0 0]); xlabel z, ylabel t, zlabel w
%         title("PPR3"); drawnow
%
% end


fprintf('\n# Table IV Data (Burgers Eq.)\n');
fprintf('# Control costs for different controllers\n');
fprintf("      Controller    &     Cost     \n")
fprintf("     %s   &  %4.3e    \n", 'Uncontrolled', costUnc)
fprintf("     %s   &  %4.3e    \n", '     LQR    ', costLQR)
fprintf("     %s   &  %4.3e    \n", '     PPR1   ', costPPR1)
fprintf("     %s   &  %4.3e    \n", '     PPR2   ', costPPR2)
fprintf("     %s   &  %4.3e    \n", '     PPR3   ', costPPR3)

