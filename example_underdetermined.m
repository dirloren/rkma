clear
rng(719)
addpath('./matlab2tikz')

writeout = true;
m = 100; % #rows
n = 500; % #columns
A = randn(m,n);
V = A.*(abs(A)>0.3); % mismatched adjoint
c = randn(m,1).*(rand(m,1)>0.0);
xhat = V'*c; % "solution" in the range of V'

%A = V + 0.4*randn(m,n); % true adjoint
b = A*xhat;

rng('shuffle')

p = sum(A.^2,2); p = p/sum(p);
scpav = diag(A*V');
pV = scpav./sum(scpav);

% check convergence condition
normv = sqrt(sum(V.^2,2));
D = diag(pV./scpav);
S = diag(normv.^2./scpav);
Z = orth(V');

rho = max(abs(eig(Z'*(eye(n) - V'*D*A)*Z)))

maxiter = 40*m;

[xV,dataV] = rkma(A,V,b,pV,maxiter);

[xA,dataA] = rkma(A,A,b,p,maxiter);

itersA = dataA.iter;
itersV = dataV.iter;

fprintf('    err       res\n');
fprintf('xV: %2.2e, %2.2e\n',norm(xhat-xV),norm(A*xV-b));
fprintf('xA: %2.2e, %2.2e\n',norm(xhat-xA),norm(A*xA-b));

errV = sqrt(sum((dataV.x - xhat*ones(1,length(itersV))).^2));
errA = sqrt(sum((dataA.x - xhat*ones(1,length(itersA))).^2));

resV = sqrt(sum((A*dataV.x - b*ones(1,length(itersV))).^2));
resA = sqrt(sum((A*dataA.x - b*ones(1,length(itersA))).^2));

clf
semilogy(itersA,errA,itersV,errV,itersV,rho.^itersV*norm(xhat))
xlabel('$k$','Interpreter','latex')
ylabel('$\|x_k-\hat x\|$','Interpreter','latex')
legend('RK','RKMA','expected rate')
if writeout
    matlab2tikz('width','\figurewidth',...
        'extraaxisoptions',['legend style={font=\scriptsize},'], ...
        '../tex/figures/example_underdetermined_error.tex');
end


%%
clf
semilogy(itersA,resA,itersV,resV)
xlabel('$k$','Interpreter','latex')
ylabel('$\|Ax_k-b\|$','Interpreter','latex')
legend('RK','RKMA')
if writeout
    matlab2tikz('width','\figurewidth',...
        'extraaxisoptions',['legend style={font=\scriptsize},'], ...
        '../tex/figures/example_underdetermined_residual.tex');
end

% check convergence condition

av = diag(A*V');
normv = sqrt(sum(V.^2,2));
D = diag(p./av);
S = diag(normv.^2./av);
Z = orth(V')'; % orthonormal basis for range of V.
M = V'*D*A + A'*D*V - A'*S*D*A;
e = eig(Z*M*Z');
lambda = min(e)
