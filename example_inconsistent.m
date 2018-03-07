clear
rng(719)
addpath('./matlab2tikz')
writeout = true;
m = 500; % #rows
n = 200; % #columns
A = randn(m,n);
V = A + randn(m,n)*0.05; % mismatched adjoint
xhat = randn(n,1);
r = 0.01*randn(m,1);
b = A*xhat + r;

%rng('shuffle')

p = sum(A.^2,2); p = p/sum(p);
scpav = diag(A*V');
pV = scpav./sum(scpav);

maxiter = 80*m;


% check convergence condition:
av = diag(A*V');
normv = sqrt(sum(V.^2,2));

D = diag(pV./av);
S = diag(normv.^2./av);

M = V'*D*A + A'*D*V - A'*S*D*A;

lambda = min(eig(M))

gamma = max(abs(r).*normv./av);

% solve with rkma
[xV,dataV] = rkma(A,V,b,p,maxiter);
itersV = dataV.iter;

resV = sqrt(sum((A*dataV.x - b*ones(1,length(itersV))).^2));
errV = sqrt(sum((dataV.x - xhat*ones(1,length(itersV))).^2));

I = eye(n);
M = I - V'*D*A;
expfinalerror = norm( (V'*D*A)\(V'*D*r) ); % is norm of sum_{l=0}^maxiter M^l*V'*D*r 

semilogy(itersV,sqrt((1-lambda/2).^(itersV)*norm(xhat)^2 + gamma^2*2/lambda),itersV,errV,itersV,expfinalerror*ones(size(itersV)))
xlabel('$k$','Interpreter','latex')
ylabel('$\|x_k-\hat x\|$','Interpreter','latex')
legend('upper bound','RKMA','expected final error')
if writeout
    matlab2tikz('width','\figurewidth',...
        'extraaxisoptions',['legend style={font=\scriptsize},'], ...
        '../tex/figures/example_inconsistent_error.tex');
end



