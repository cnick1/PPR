function [f, g, h, x0] = getSystem31(reduced, initialCondition)
%getSystem31 Generates a finite difference model for a 1D Fokker-Planck equation from [1]
%
%   Usage:  [f,g,h] = getSystem31()
%
%   Inputs:    reduced - boolean; whether to use the n=99 FOM (false) or n=25 ROM (true)
%     initialCondition - '1' or '2', for test case 1 or 2 from [1]
%
%   Outputs:     f,g,h - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%
%   Description: The controlled Fokker-Planck equation, as in [1], is given by
%
%     ∂ρ/∂t = νΔρ + ∇·(ρ∇G) + Σᵐⱼ₌₁ u_j ∇·(ρ∇α_j) in Ω×(0,∞),
%     0 = (ν∇ρ + ρ∇G)·n on Γ×(0,∞),
%     ρ(x,0) = ρ0(x) in Ω.
%
%   For more details, see [1]. The finite difference approximation returned 
%   here is based on code provided by the authors of [1] via email
%   correspondence. All of this code was written by them, and is just
%   repackaged here for convenient usage with the standardized PPR package
%   syntax.
%
%   References: [1] T. Breiten, K. Kunisch, and L. Pfeiffer, “Numerical study of
%               polynomial feedback laws for a bilinear control problem,”
%               Mathematical Control & Related Fields, vol. 8, no. 3, pp.
%               557–582, 2018, doi: 10.3934/mcrf.2018023.
%
%%
if nargin < 2
    initialCondition = '2';
    if nargin < 1
        reduced = false;
    end
end
rng(1)

nx = 100;
xa = -6;
xb = 6;
nu = 1;
hx = (xb-xa)/(nx+1);
xd = xa+hx:hx:xb-hx;

% backward Kolmogorov equation
% d/dt y = nu*Delta*y - < nabla y , nabla W >
% < nabla y , normal > = 0
%
W = @(x) ( ((0.5*x.^2 - 15).*x.^2 + 119).*x.^2 + 28*x + 50) / 200;
dWdx  = @(x) (3*x.^5)/200 - (3*x.^3)/10 + (119*x)/100 + 7/50;

ex = ones(nx,1);
I = speye(nx);
Laplace_x = 1/hx^2*spdiags([ex -2*ex ex], -1:1, nx, nx);
Grad_x_fw = 1/hx*spdiags([-ex ex],0:1,nx,nx);
Grad_x_bw = 1/hx*spdiags([-ex ex],-1:0,nx,nx);
Grad_x = 1/(2*hx)*spdiags([-ex 0*ex ex],-1:1,nx,nx);

dWx = spdiags(reshape(dWdx(xd)',nx,1),0,nx,nx);

e1x = I(:,1);
enx = I(:,nx);

O = sparse(nx,nx);


A = nu*Laplace_x- min(dWx,O)*Grad_x_fw - max(dWx,O)*Grad_x_bw;
A = A + nu/hx^2*(e1x*e1x') + nu/hx^2*(enx*enx') - 1/hx*min(dWx,0)*(enx*enx')...
    + 1/hx*max(dWx,0)*(e1x*e1x') ;

e = hx*ones(nx,1);

A = A';
n = size(A,1);%
f = @(x,gap) (xa+gap < x && x < xa+2*gap)*1/(xb-xa)*(2/gap^2*(x-xa-gap)^3 ...
    -5/gap^3*(x-xa-gap)^3*(x-xa-2*gap)+9/gap^4*(x-xa-gap)^3*(x-xa-2*gap)^2) ...
    + (xa+2*gap <= x && x <= xb-2*gap)*1/(xb-xa)*(x-xa) ...
    + (xb-2*gap < x && x < xb -gap)*(1-2*gap/(xb-xa)+1/(xb-xa)*(x-xb+2*gap) ...
    +1/(gap^2*(xb-xa))*(x-xb+2*gap)^3-4/(gap^3*(xb-xa))*(x-xb+2*gap)^3*(x-xb+gap)...
    + 9/(gap^4*(xb-xa))*(x-xb+2*gap)^3*(x-xb+gap)^2)+ (xb-gap <= x && x<= xb);

%%
alpha = zeros(nx,1);
for i = 1:nx
    alpha(i) = f(xa+(i-1)*hx,1e-1)-0.5;
end

%%

Diffx = Grad_x - 1/(2*hx)*(e1x*e1x') + 1/(2*hx)*(enx*enx');
% Computing control matrices N,B
Nt = -spdiags(Diffx*alpha,0,nx,nx)*Grad_x;
Nt = Nt + 1/(2*hx)*spdiags(Diffx*alpha,0,nx,nx)*(e1x*e1x') - ...
    1/(2*hx)*spdiags(Diffx*alpha,0,nx,nx)*(enx*enx');
N = Nt';

pow_err = inf;
rho = randn(n,1);
while pow_err > 1e-10
    rho_old = rho;
    rho_new = (A-speye(n))\rho;
    rho = rho_new/(e'*rho_new);
    pow_err = norm(rho-rho_old)/norm(rho_old);
end
B = N*rho;
orig_bil = @(t,x,A) A*x;

%% projection onto stable subspace
m = size(A,1)-1;
n = size(A,1);
R = [speye(m),rho(1:m);-ones(1,m),rho(n)];
iR = inv(R);
% iR = [speye(m),sparse(m,1);e(1:n)']-[rho(1:m);0]*e';
Q = [speye(m);zeros(1,m)];
C = speye(nx);

rA = sparse(round(1e12*(iR*A*R))/1e12);
rN = sparse(round(1e12*(iR*N*R))/1e12);
rB = round(1e12*(iR*B))/1e12;
rC = C*R*Q*sqrt(hx);

tA = rA(1:m,1:m);
tN = rN(1:m,1:m);
tB = rB(1:m,:);
tC = rC(:,1:m);

%% Construct initial condition
switch initialCondition
    case '1'
        rng(1)
        x0 = rand(nx,1);
        x0 = x0/(e'*x0);
        x0 = x0+rho;
        x0 = x0/(e'*x0);
    case '2'
        x0 = 1*exp(-(xd').^2/(2*0.5^2));
        x0 = x0/(e'*x0);
end

%% Construct f,g,h cell arrays
% [1] uses these "transformed" operators
if reduced % Return ROM
    [ Ar,Nr,Br,Cr,Vr,Wr ] = compute_ROM( tA,tN,tB,tC,rel_tol_lyap,trunc_tol,red_tol );
    f = {Ar}; g = {Br,Nr}; h = {Cr};
    x0 = Wr'*(Q'*(R\(x0-rho)));
else % Return FOM
    f = {tA}; g = {tB,tN}; h = {tC};
    x0 = Q'*(R\(x0-rho)); % Transform IC to model coordinates
end

end

function [ Ar,Nr,Br,Cr,Vr,Wr ] = compute_ROM( A,N,B,C,rel_tol_lyap,trunc_tol,red_tol )
%COMPUTE_ROM Summary of this function goes here
%   Detailed explanation goes here
% Note: this function was written by the authors of [1], and is reproduced
% here *as-is* without modification.

P = lyap(A,B*B');
res = norm(A*P+P*A'+N*P*N'+B*B','fro')/norm(B*B','fro')
while res > rel_tol_lyap
    P = lyap(A,B*B'+N*P*N');
    res = norm(A*P+P*A'+N*P*N'+B*B','fro')/norm(B*B','fro')
end

Q = lyap(A',C'*C);
res = norm(A'*Q+Q*A+N'*Q*N+C'*C,'fro')/norm(C'*C,'fro')
while res > rel_tol_lyap
    Q = lyap(A',C'*C+N'*Q*N);
    res = norm(A'*Q+Q*A+N'*Q*N+C'*C,'fro')/norm(C'*C,'fro')
end

P = 0.5*(P+P');
Q = 0.5*(Q+Q');
[UP,DP] = eig(P);
[UQ,DQ] = eig(Q);
Ptemp.L = UP;
Ptemp.D = DP;
Qtemp.L = UQ;
Qtemp.D = DQ;
X = truncation(Ptemp,trunc_tol);
Y = truncation(Qtemp,trunc_tol);
S = (X.L*sqrt(X.D))';
RR = (Y.L*sqrt(Y.D))';
[U,Sigma,V] = svd(S*RR',0); 
% figure(1)
% semilogy(diag(Sigma)/Sigma(1,1),'o')
% r = input('Specify dimension of ROM: ');
% close all
r = find(diag(Sigma)/Sigma(1,1) < red_tol,1,'first');
% r = 8;
U = U(:,1:r);
Sigma = Sigma(1:r,1:r);


V = V(:,1:r);
%
Vr = S'*U*Sigma^-0.5;
Wr = RR'*V*Sigma^-0.5;
%
Ar = Wr'*A*Vr;
Nr = Wr'*N*Vr;
Br = Wr'*B;
Cr = C*Vr;

end

