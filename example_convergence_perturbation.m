rng(719)
rng('shuffle')
addpath('./matlab2tikz')
writeout = false;
m = 500; % #rows
n = 200; % #columns
A = randn(m,n);
V = A.*(abs(A)>0.5); % mismatched adjoint
xhat = randn(n,1);
b = A*xhat;

p = sum(A.^2,2); p = p/sum(p);
scpav = diag(A*V');
pV = scpav./sum(scpav);

maxiter = 40*m;


% check convergence condition:
av = diag(A*V');
normv = sqrt(sum(V.^2,2));

D = diag(pV./av);
S = diag(normv.^2./av);

M = V'*D*A + A'*D*V - A'*S*D*A;

lambda = min(eig(M));
1-lambda

M = eye(n) - V'*D*A;
rho = max(abs(eig(M)))

norm(M)

%return

% solve with rkma
[xV,dataV] = rkma(A,V,b,pV,maxiter);
[xA,dataA] = rkma(A,A,b,p,maxiter);
itersV = dataV.iter;
itersA = dataA.iter;

errA = sqrt(sum((dataA.x - xhat*ones(1,length(itersA))).^2));
resA = sqrt(sum((A*dataA.x - b*ones(1,length(itersA))).^2));
errV = sqrt(sum((dataV.x - xhat*ones(1,length(itersV))).^2));
resV = sqrt(sum((A*dataV.x - b*ones(1,length(itersV))).^2));


semilogy(itersA,errA,itersV,errV,itersV,rho.^itersV*norm(xhat))
xlabel('$k$','Interpreter','latex')
ylabel('$\|x_k-\hat x\|$','Interpreter','latex')
legend('RK','RKMA','expected rate')
if writeout
    matlab2tikz('width','\figurewidth',...
        'extraaxisoptions',['legend style={font=\scriptsize},'], ...
                        '../tex/figures/example_convergence_perturbation_error.tex');
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
                        '../tex/figures/example_convergence_perturbation_residual.tex');
end

lambdaA = min(eig(2*A'*D*A - A'*S*D*A))
