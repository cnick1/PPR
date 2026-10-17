% function runExample31_scalability()
%runExample31_scalability Runs the 1D Fokker-Planck control example to compare with [1].
%
%   Usage:  runExample31_scalability()
%
%   Description: The controlled Fokker-Planck equation, as in [1], is given by
%
%     ∂ρ/∂t = νΔρ + ∇·(ρ∇G) + Σᵐⱼ₌₁ u_j ∇·(ρ∇α_j) in Ω×(0,∞),
%     0 = (ν∇ρ + ρ∇G)·n on Γ×(0,∞),
%     ρ(x,0) = ρ0(x) in Ω.
%
%   For more details on the model, see [1]. The finite difference model
%   approximation used here is based on code provided by the authors of
%   [1] via email correspondence.
%
%   Upon discretization and the coordinate changes necessary to center the
%   model, we obtain a finite-dimensional bilinear model of the form
%
%           ẋ = A x + B u + G₁ (x ⊗ u)
%
%   for which we wish to solve an optimal control problem minimizing the
%   cost function
%
%           J(x,u) = ½∫ xᵀ Q x + uᵀ R u dt
%
%   The solution is given by the value function V(x) which solves the HJB
%   PDE, and the optimal control is then
%
%           u(x) = -R⁻¹ gᵀ(x) 𝜕V(x)/𝜕x
%
%   The value function approximation in [1] is in principle identical to
%   that computed in the PPR package, i.e. a Taylor approximation computed
%   via Al'brekht's method. The details of their implementation are
%   slightly different, most notably their use of a quadrature-based k-way
%   Lyapunov solver as opposed to the Bartels-Stewart type solvers the PPR
%   package uses.
%
%   However, the optimal control approximation proposed in [1] is slightly
%   different than the one proposed in the PPR package. Due to the
%   nonlinear input map g(x), it is possible to compute a polynomial
%   feedback law with *some* higher-order terms than the value function
%   approximation. For example, a quadratic value function approximation
%   with Riccati solution V₂ defines a linear feedback law in a Taylor
%   sense: u₁(x) = -R⁻¹ Bᵀ V₂ x. A quadratic *Taylor* approximation of the
%   feedback law would require a cubic value function approximation, as
%   some---but not all---of the quadratic terms in the feedback law contain
%   the cubic value function coefficient. However, due to the nonlinear
%   input g(x), some of the quadratic feedback law terms rely only on lower
%   value function coefficients. One can, therefore, construct *A*
%   quadratic polynomial feedback law using only V₂, simply by including
%   the entire polynomial input map:
%
%       u₂(x) = -R⁻¹ gᵀ(x) V₂ x
%
%   This is what [1] computes, and it is NOT a Taylor approximation (nor do
%   they explicitly claim it is). The PPR package can compute these
%   approximations as well; for bilinear systems, it is straightforward to
%   adapt the code to compute one additional K coefficient. However, it is
%   more efficient to compute g(x), vᵈᵀ(x ⊗ ... ⊗ x), and then multiply
%   them together, so we do this here instead.
%
%   In this script, we will compare scalability of the two approaches.
%
%   References: [1] T. Breiten, K. Kunisch, and L. Pfeiffer, “Numerical study of
%               polynomial feedback laws for a bilinear control problem,”
%               Mathematical Control & Related Fields, vol. 8, no. 3, pp.
%               557–582, 2018, doi: 10.3934/mcrf.2018023.
%
%   Part of the PPR repository.
%%
fprintf('Running Example 31\n')

%% Build Table I
%  Computational performance of the energy function approximations.
%  Since the initial times are so short, we average nTest times
degree = 3; nTest = 10;
nd = []; times = []; times_breiten=[]; vOfX = []; 
nxs = [16, 32, 64];%, 128, 256, 512, 1024];
for nx = nxs
    if nx == 128; nTest = 1; end % only average over the ones with very fast runtime
    [f, g, h, x0] = getSystem31(false, '2',nx); x0 = x0./2;
    q = h{1}.'*h{1}; R = 1e-4;

    tic; for i = 1:nTest, v = breitenValueFunction(f{1},g{1},g{2},q,R,50,degree); end, tt_breiten = toc / nTest;
    tic; for i = 1:nTest, [v, ~] = ppr(f, g, q, R, degree); end, tt = toc / nTest;

    nd = [nd, length(v{end})]; times = [times, tt]; times_breiten = [times_breiten, tt_breiten];
    vOfX = [vOfX, 0.5 * kronPolyEval(v, x0)];
end


fprintf('Writing data to plots/example31_convergenceData_d3.dat \n')
fileID = fopen('plots/example31_convergenceData_d3.dat', 'w');
fprintf(fileID, '# Table I Data\n# Fokker-Planck model, convergence and scalability results; d=%d \n', degree);
fprintf(fileID, '    n &        n^%d &   CPU-sec (PPR)   &   CPU-sec (Breiten)   &      V(x0)      \n', degree);
for i = 1:length(nxs)
    fprintf(fileID, '%5d & %10d &    %8.5e    &    %8.5e    &    %8.5e       \n', nxs(i)-1, nd(i), times(i), times_breiten(i), vOfX(i));
end
fclose(fileID);

type('plots/example31_convergenceData_d3.dat');


%%
degree = 4; nTest = 10;
nd = []; times = []; times_breiten=[]; vOfX = []; 
nxs = [16, 32];%, 64, 128];
for nx = nxs
    if nx == 64; nTest = 1; end % only average over the ones with very fast runtime
    [f, g, h, x0] = getSystem31(false, '2',nx); x0 = x0./2;
    q = h{1}.'*h{1}; R = 1e-4;

    tic; for i = 1:nTest, v = breitenValueFunction(f{1},g{1},g{2},q,R,50,degree); end, tt_breiten = toc / nTest;
    tic; for i = 1:nTest, [v, ~] = ppr(f, g, q, R, degree); end, tt = toc / nTest;

    nd = [nd, length(v{end})]; times = [times, tt]; times_breiten = [times_breiten, tt_breiten];
    vOfX = [vOfX, 0.5 * kronPolyEval(v, x0)];
end


fprintf('Writing data to plots/example31_convergenceData_d4.dat \n')
fileID = fopen('plots/example31_convergenceData_d4.dat', 'w');
fprintf(fileID, '# Table I Data\n# Fokker-Planck model, convergence and scalability results; d=%d \n', degree);
fprintf(fileID, '    n &        n^%d &   CPU-sec (PPR)   &   CPU-sec (Breiten)   &      V(x0)      \n', degree);
for i = 1:length(nxs)
    fprintf(fileID, '%5d & %10d &    %8.5e    &    %8.5e    &    %8.5e        \n', nxs(i)-1, nd(i), times(i), times_breiten(i), vOfX(i));
end
fclose(fileID);

type('plots/example31_convergenceData_d4.dat');


%% Reduced
% degree = 3; nTest = 10;
% nd = []; times = []; times_reduced=[]; vOfX = []; vrOfX = []; 
% nxs = [16, 32, 64, 128, 256, 512];
% rxs = [15,15,15,15,15,15,15];
% j=0;
% for nx = nxs
%     j=j+1;
%     if nx == 64; nTest = 1; end % only average over the ones with very fast runtime
%     [f, g, h, x0] = getSystem31(false, '2',nx); x0 = x0./2;
%     q = h{1}.'*h{1}; R = 1e-4; options.reducedDimension = rxs(j);
% 
%     tic; for i = 1:nTest, [vr, ~] = ppr(f, g, q, R, degree, options); end, tt_r = toc / nTest;
%     tic; for i = 1:nTest, [v, ~] = ppr(f, g, q, R, degree); end, tt = toc / nTest;
% 
%     nd = [nd, length(v{end})]; times = [times, tt]; times_reduced = [times_reduced, tt_r];
%     vOfX = [vOfX, 0.5 * kronPolyEval(v, x0)];
%     vrOfX = [vrOfX, 0.5 * kronPolyEval(vr, x0)];
% end
% 
% 
% fprintf('Writing data to plots/example31_convergenceData_d3_r.dat \n')
% fileID = fopen('plots/example31_convergenceData_d3_r.dat', 'w');
% fprintf(fileID, '# Table I Data\n# Fokker-Planck model, convergence and scalability results; d=%d \n', degree);
% fprintf(fileID, '    n &        n^%d &   CPU-sec (PPR)   &   CPU-sec (PPR Red.)   &      V(x0)    &      Vr(x0)      \n', degree);
% for i = 1:length(nxs)
%     fprintf(fileID, '%5d & %10d &    %8.5e    &    %8.5e    &    %8.5e    &    %8.5e        \n', nxs(i)-1, nd(i), times(i), times_reduced(i), vOfX(i), vrOfX(i));
% end
% fclose(fileID);
% type('plots/example31_convergenceData_d3_r.dat');
% 
% %% Reduced
% degree = 4; 
% nd = []; times_reduced=[]; vrOfX = []; 
% nxs = [16, 32, 64, 128];
% rxs = [15,15,15,15,15];
% j=0;
% for nx = nxs
%     j=j+1;
%     [f, g, h, x0] = getSystem31(false, '2',nx); x0 = x0./2;
%     q = h{1}.'*h{1}; R = 1e-4; options.reducedDimension = rxs(j);
% 
%     tic; [vr, ~] = ppr(f, g, q, R, degree, options); tt_r = toc;
% 
%     nd = [nd, (nx-1)^degree]; times_reduced = [times_reduced, tt_r];
%     vrOfX = [vrOfX, 0.5 * kronPolyEval(vr, x0)];
% end
% 
% fprintf('Writing data to plots/example31_convergenceData_d4_r.dat \n')
% fileID = fopen('plots/example31_convergenceData_d4_r.dat', 'w');
% fprintf(fileID, '# Table I Data\n# Fokker-Planck model, convergence and scalability results; d=%d \n', degree);
% fprintf(fileID, '    n &        n^%d &   CPU-sec (PPR)   &   CPU-sec (PPR Red.)   &      V(x0)    &      Vr(x0)      \n', degree);
% for i = 1:length(nxs)
%     fprintf(fileID, '%5d & %10d &    %8.5e    &    %8.5e    &    %8.5e    &    %8.5e        \n', nxs(i)-1, nd(i), times(i), times_reduced(i), vOfX(i), vrOfX(i));
% end
% fclose(fileID);
% type('plots/example31_convergenceData_d4_r.dat');
% 
% %%
% type('plots/example31_convergenceData_d4.dat');
% fprintf('    n &        n^%d &   CPU-sec (PPR reduced)   &           &     Vr(x0)      \n', degree);
% for i = 1:length(nxs)
%     fprintf('%5d & %10d &    %8.5e    &                   &    %8.5e        \n', nxs(i)-1, nd(i), times_reduced(i), vrOfX(i));
% end




function y = gradTermEval(v,x,k)
% Compute the k-1th degree term of the gradient of the value function (the
% derivative term corresponding to the kth degree term in the value
% function). Due to the symmetry of vk, we can do this by reshaping and
% multiplying by x, k-1 times, contracting the tensor Vk down to a vector
% as we go.
n = size(x, 1);
y = v{k};
for i = k:-1:2
    % Using kron-vec identity (level-3 vs level-2 BLAS)
    y = reshape(y,n^(i-1),n) * x;
end
end

%% The code below was provided by Breiten et al.
function gain = computeGain(T,x)
% Note: this is moreso just the gradient of the value function, not really
% the gain matrix. This is their implementation, which is not as efficient
% as what we do (you can recycle xt instead of starting from scratch every
% time, and instead of forming x^k you can just multiply with x k-times
% (which is what gradTermEval above does).
% Their T is akin to our v, but our v's are symmetric vectors (vectorized
% symmetric tensors) whereas it is not clear without digging further what
% properties their T's satisfy, so I will leave their implementation intact.
n = size(x, 1);
gain = 0;
for k = 2:size(T,2)
    xt = x';
    for j = 1:k-2
        xt = kron(x',xt);
    end
    gain = gain + 1/factorial(k-1)*xt*reshape(T{k},n^(k-1),n);
end
end

function T = breitenValueFunction(A,B,G1,q,R,quad_no,degree)
% This is a cleaned up version of the code from [1], packaged as a function
T{2} = icare(full(A),full(B),full(q),R);
Api = A-R\B*B'*T{2};

for j = 3:degree
    rhs  = set_up_rhs(T,B,G1,j,R);

    % Approximation by Quadrature
    T{j} = solv_kron_ls( Api',rhs,j,quad_no );
end

end

function [ rhs ] = set_up_rhs( T,B,N,k,gamma )
%SET_UP_RHS Summary of this function goes here
%   Detailed explanation goes here
n = size(B,1);
u{1} = mat_mult( T{2},2,n,1,B' );
for j = 2:k-2
    u{j} = reshape(mat_mult( T{j+1},j+1,n,1,B' ),n*ones(1,j));
    for l = 1:j
        u{j} = u{j} + mat_mult( T{j},j,n,l,N' );
    end
end

u{k-1} = 0;
for l = 1:k-1
    u{k-1} = u{k-1} + mat_mult( T{k-1},k-1,n,l,N' );
end

rhs = 0;
for i = 1:k-1
    v = kron(u{i}(:),u{k-i}(:));
    v = reshape(v,n*ones(1,k));

    p = nchoosek(1:k,i);
    clear q
    for l = 1:size(p,1)
        q(l,:) = setdiff(1:k,p(l,:));
    end
    PQ = [q,p];
    for j = 1:size(p,1)
        perv = ipermute(v,PQ(j,:));
        rhs = rhs + perv;
    end
end
rhs = rhs/(2*gamma);

end

function [ x ] = solv_kron_ls( A,f,d,k )
%COMPUTE_INVERSE Summary of this function goes here
%   Detailed explanation goes here

lam_min_A = -max(real(eig(A)));
lam_min = d*lam_min_A;
x = 0;
hst = pi/sqrt(k);
for j = -k:k
    tj = log(exp(j*hst)+sqrt(1+exp(2*j*hst)));
    wj = hst/sqrt(1+exp(-2*j*hst));
    D = expm(2*tj/lam_min*A);
    x = x-2*wj/lam_min*mat_mult_all( f,D,d );
    % j
end


end

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

while (~test)&&(k<1000)
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

