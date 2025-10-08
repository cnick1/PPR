function runExample26(degree, transformed)
%runExample26 Runs 2D inverted pendulum example (normalized). 
%   Plots a) the value function, b) the HJB residual, and c) closed-loop phase portraits.
%
%   Usage:  runExample26()
%
%   Inputs:    degree - degree of value function approximation
%         transformed - whether to use model (1) (false) or (2) (true)
%        
%   Description: The basic equation for the inverted pendulum is
%
%               ẍ - sin(x) = u(t)
%
%   This can be put in first-order form with x₁ = x, x₂ = ẋ as   
%
%             ẋ₁ = x₂
%      (1)    ẋ₂ = sin(x₁) + u(t)
%              y = x₁
%
%   Alternatively, we can transform with the following (bijective)
%   coordinate transformation
%
%              x = Φ(z)   = [z₁; z₂ - 2 sin(z₁/2)]
%              z = Φ⁻¹(x) = [x₁; x₂ + 2 sin(x₁/2)]
%
%   to obtain the transformed dynamical system
%
%             ̇z₁ = z₂ - 2 sin(z₁/2)
%      (2)    ̇z₂ = z₂ cos(z₁/2) + u(t)
%              y = z₁
% 
%   Now, we can compute a controller u(x) = K(x) using PPR. Alternatively, 
%   we can compute a controller u(z) = K(z) = K(Φ⁻¹(x)). Is there a
%   difference? Is one better than the other? 
%
%   Hypothesis: Computing the controller in the z coordinates will be
%   easier with polynomials than in the x coordinates, potentially allowing
%   global stabilization. In other words, whereas u(x) is only locally 
%   stabilizing with PPR, I am hoping that u(z) is globally stabilizing. 
% 
%   References: [1] 
%
%   Part of the PPR repository.
%%
if nargin < 2 
    transformed = false;
    if nargin < 1 
    degree = 6;
    end
end

fprintf('Running Example 26\n')

reducedOrder = 1;
[f, g, ~] = getSystem26(degree-1,transformed);

F = @(x) [x(2); sin(x(1))];
G = @(x) g{1};

%% Open-loop phase portrait
% plotPhasePortrait(1, @(x) 0, true)
% plotPhasePortrait(1, @(x) 0, false)

%% Value function & HJB Residual plots
%  Compute the polynomial approximations to the value function
q = 0; r = 1; tic; 
[v, K, ~] = ppr(f, g, q, r, degree);
fprintf("Completed ppr() in %2.2f seconds. \n", toc)

nX = 301; nY = nX; xLim = pi; yLim = 5;
xPlot = linspace(-xLim, xLim, nX); yPlot = linspace(-yLim, yLim, nY); [X, Y] = meshgrid(xPlot, yPlot);

HJBResidual = zeros(nY, nX); valueFunction = zeros(nY, nX); 

uPPR = @(x) kronPolyEval(K, x);
V = @(x) 0.5 * kronPolyEval(v,x);
dVdx = @(x) 0.5 * kronPolyDerivEval(v, x);

for i = 1:nY
    for j = 1:nX
        x = [X(i, j); Y(i, j)];
        valueFunction(i, j) = V(x);
        HJBResidual(i, j) = dVdx(x)*F(x) - 1/2*dVdx(x)*G(x)*G(x).'*dVdx(x).';
        if valueFunction(i, j) < 0
            valueFunction(i, j) = NaN;
        end
    end
end

fig1 = figure;
contourf(X, Y, valueFunction, 16, 'w'); hold on;
xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16); xticks([-pi, 0, pi]); xticklabels({'-\pi', '0', '\pi'})
load(fullfile('utils', 'YlGnBuRescaled.mat')); colormap(flip(YlGnBuRescaled))
% clim([0 4e4])
drawnow

axis off
fprintf('Exporting figure to: \n     plots/example26_valueFun_d%i_polynomial%i.pdf\n', degree, degree-1)
exportgraphics(fig1, sprintf('plots/example26_valueFun_d%i_polynomial%i.pdf', degree, degree-1), 'ContentType', 'vector', 'BackgroundColor', 'none');
colorbar

fig2 = figure;
pcolor(X, Y, log10(abs(HJBResidual))); shading interp;
xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex');
set(gca, 'FontSize', 16); xticks([-pi, 0, pi]); xticklabels({'-\pi', '0', '\pi'})
load('utils\YlGnBuRescaled.mat'); colormap(flip(YlGnBuRescaled))
clim([-3 9])
drawnow

axis off
fprintf('Exporting figure to: \n     plots/example26_valueFun-HJB-Error_d%i_polynomial%i.pdf\n', degree, degree-1)
exportgraphics(fig2, sprintf('plots/example26_valueFun-HJB-Error_d%i_polynomial%i.pdf', degree, degree-1), 'ContentType', 'vector', 'BackgroundColor', 'none');


%% Closed-loop phase portraits
plotPhasePortrait(degree, uPPR, transformed)

%% Estimate region of attraction with Zubov 
% uLQR = @(x) kronPolyEval(K, x, 1);
% 
% cmap = [1, 0.6, 0.6;  % Light red (V < 0)
%         0.6, 1, 0.6]; % Light green (V > 0)
% 
% Vlin = zeros(nY, nX); 
% Vdotlin = zeros(nY, nX); 
% 
% for i = 1:nY
%     for j = 1:nX
%         x = [X(i, j); Y(i, j)];
% 
%         Vlin(i,j) = kronPolyEval(v,x);
%         Vdotlin(i,j) = kronPolyDerivEval(v,x) * (F(x) + G(x)*uPPR(x));
%         % Vdotlin(i,j) = kronPolyDerivEval(v,x, 2) * (F(x) + G(x)*uLQR(x));
% 
%     end
% end
% 
% figure; contourf(X, Y, Vlin); colorbar; title("Lyapunov function V(x)")
% figure; contourf(X, Y, Vdotlin); colorbar; title("Lyapunov function Lie derivative dV(x)/dt")
% figure; pcolor(X, Y, Vdotlin); shading interp; colorbar; title("Lyapunov function Lie derivative dV(x)/dt")
% caxis([-1000 100])
% % Plot region of attraction using nonlinear Lyap estimate 
% figure; hold on; xlabel('$x_1$', 'interpreter', 'latex'); ylabel('$x_2$', 'interpreter', 'latex'); title('R_A estimate using polynomial Lyapunov function')
% % Define colors: light green for V > 0, red for V < 0
% colormap(cmap)
% 
% % Plot region where V is positive definite 
% [c,h] = contourf(X, Y, Vlin,[-inf 0]); h.FaceAlpha = 0.25;
% 
% % Plot region where Vdot is negative definite 
% [c,h] = contourf(X, Y, -Vdotlin,[-inf 0]); h.FaceAlpha = 0.25;
% 
% % Plot largest Lyapunov level set that fits in the region 
% ax1 = gca; xlim([-xLim xLim]); ylim([-yLim yLim]);
% ax2 = axes; xlim([-xLim xLim]); ylim([-yLim yLim]);
% 
% clin = 0.997;
% contour(X, Y, 1-exp(-Vlin), clin*[1 1],'k'); ax2.Visible = 'off';    



end



function plotPhasePortrait(degree, u, transformed)

if nargin < 3
    transformed = false;
end

% Write dynamics to file
fileID = fopen('examples\example26ExportedControlFunction.py', 'w'); 

% if transformed
%     % Write z = Φ⁻¹(x) = [x₁; x₂ + 2 sin(x₁/2)]
%     % PhiInv = @(x) [x(1); x(2) + 2 * sin(x(1)/2)];
%     PhiInv = 'PhiInv2(x1, x2)';
% else
%     % Write z = Φ⁻¹(x) = [x₁; x₂]
%     % PhiInv = @(x) [x(1); x(2)];
%     PhiInv = 'PhiInv1(x1, x2)';
% end

% Define control law function declaration
funcStr = 'def u(z1, z2):\n    ';

% Now add the z = Φ⁻¹(x) transformation 
% syms x1 x2
% funcStr = [funcStr '[z1, z2] = ' PhiInv '\n    '];

% Now add the return statement where you print the actual control law in z
syms z1 z2
funcStr = [funcStr 'return ' char(vpa(u([z1;z2])))];

% Replace all the matlab stuff with python stuff
funcStr = strrep(funcStr, '^', '**'); % Replace MATLAB's ^ operator with Python's ** operator
funcStr = strrep(funcStr, 'sin', 'np.sin'); % Replace MATLAB's sin with np.sin
funcStr = strrep(funcStr, ';', ','); % Replace MATLAB's ; with , for the array

fprintf(fileID, funcStr);

fclose(fileID);

pyenv('ExecutionMode', 'OutOfProcess'); % Optional: You can switch to 'InProcess' mode if needed
if transformed
    pyrunfile("examples\plotExample26_transformed.py", controlDegree=degree-1)
    open('plots/example26_phasePortrait.pdf')
    pyrunfile("examples\plotExample26_retransformed.py", controlDegree=degree-1)
    open('plots/example26_phasePortrait.pdf')

else
    pyrunfile("examples\plotExample26.py", controlDegree=degree-1)
    open('plots/example26_phasePortrait.pdf')
end
terminate(pyenv);
end

