function runExample28()
%runExample28 Runs the 2D linear heat equation FEM control example.
%
%   Usage:  runExample28()
%
%   Description: 
%
%
%   Now, we can compute a controller u(x) = K(x) using PPR.
%
%   References: [1]
%
%   Part of the PPR repository.
%%
fprintf('Running Example 28\n')

% Get dynamics
numElements = 16; nx = numElements+1; ny = nx;
n = nx*ny; m = 4;
[f, g, h, xyg] = getSystem28(numElements);
F = @(x) f{1}*x; G = @(x) g{1};

% Get value function/controller
q = 1; R = 1; degree = 2;
[~, K] = ppr(f, g, q, R, degree);

uPPR = @(x) kronPolyEval(K, x, degree=degree-1);
% uPPR = @(x) [0;0;1;0]; 

%% Simulate closed-loop system
X = reshape(xyg(:,1),nx,ny); 
Y = reshape(xyg(:,2),nx,ny);
x0 = sin(4*pi*X) + cos(3*pi*Y)+2;
x0 = x0(:); tspan = [0 5000];
    
[t, U] = ode45(@(t, x) F(x) + G(x)*uPPR(x), tspan, x0);

%% Plot solution
figure('Position', [311.6667 239.6667 1.0693e+03 573.3333]); 
for i=1:100:length(t)
    Z = reshape(U(i,:),nx,ny);

    subplot(1,2,1); grid on;
    [c,h] = contour(X,Y,Z); clabel(c,h)
    xlabel('x, m'); ylabel('y, m'); axis equal

    subplot(1,2,2); surfc(X,Y,Z); zlim([-1 5])
    xlabel('x, m'); ylabel('y, m');

    drawnow
end

end



