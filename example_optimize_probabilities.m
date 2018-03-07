clear
rng(178392)
addpath('./matlab2tikz')
writeout = false;
% Build a two matrices
n = 100;
m = 300;

A = randn(m,n);
mask = (rand(m,n)>0.05);%(abs(A)>0.1);
A = diag(2./(sqrt(1:m)+5))*A;
V = A.*mask;

%V = A.*(abs(A)>0.0);
%V = A;

%%
norma = sqrt(sum(A.^2,2));
normv = sqrt(sum(V.^2,2));
scpav = sum(A.*V,2);
S = diag(normv.^2./scpav);


%% optimize p by supergradient ascent on lambda
% initialize with uniform probability
p = ones(m,1)/m;
% initialize with rownorms
%p = scpav/sum(scpav);
cf = [];
for k=1:600
  % form matrix 
  M = V'*diag(p./scpav)*A + A'*diag(p./scpav)*V - A'*S*diag(p./scpav)*A;
  % compute eigendecomp
  [Q,Lambda] = eig(M);
  % find smallest eigenvalue and vector
  [lambda,j] = min(diag(Lambda));
  % visualize
  cf(k) = lambda;
  subplot(1,2,1)
  plot(p);
  %title(num2str(lambda))
  subplot(1,2,2)
  plot(1:k,cf);
  drawnow
  %pause
  % do subgradient ascent
  x = Q(:,j);
  %t = 1/(mod(k,50)^1+10);
  %t = 1/(floor(k/1) + 10);
  t = 1/(k+5);
  p = p + t*(((2*V-S*A)*x).*(A*x))./scpav;
  % project back to simplex
  p = projsplx(p);
  %p = p/sum(p);
end
poptl = p;

%% optimize p by subgradient descent on spectral norm of V^TDA
% initialize with uniform probability
p = ones(m,1)/m;
% initialize with rownorms
%p = scpav/sum(scpav);
cf = [];
I = eye(n);
for k=1:200
  % form matrix 
  M = I - V'*diag(p./scpav)*A;
  % compute dominating singular value and vectors (assuming it is unique)
  [u,sigma,v] = svds(M,1);
  g = -(A*u).*(V*v)./scpav;
  % visualize
  cf(k) = sigma;
  subplot(1,2,1)
  plot(p);
  %title(num2str(lambda))
  subplot(1,2,2)
  plot(1:k,cf);
  drawnow
  %pause
  % do subgradient ascent
  t = 1/(k+5);
  p = p - t*g;
  % project back to simplex
  p = projsplx(p);
  %p = p/sum(p);
end
poptn = p;

%% solve with different probabilities

%right hand side
xhat = randn(n,1);
b = A*xhat;

punif = ones(m,1)/m;
prow = scpav/sum(scpav);

maxiter = 80*n;
[xoptl,dataoptl] = rkma(A,V,b,poptl,maxiter);
[xoptn,dataoptn] = rkma(A,V,b,poptn,maxiter);
[xunif,dataunif] = rkma(A,V,b,punif,maxiter);
[xrow,datarow] = rkma(A,V,b,prow,maxiter);

itersoptl = dataoptl.iter;
erroptl = sqrt(sum((dataoptl.x - xhat*ones(1,length(itersoptl))).^2));
resoptl = sqrt(sum((A*dataoptl.x - b*ones(1,length(itersoptl))).^2));

itersoptn = dataoptn.iter;
erroptn = sqrt(sum((dataoptn.x - xhat*ones(1,length(itersoptn))).^2));
resoptn = sqrt(sum((A*dataoptn.x - b*ones(1,length(itersoptn))).^2));

itersunif = dataunif.iter;
errunif = sqrt(sum((dataunif.x - xhat*ones(1,length(itersunif))).^2));
resunif = sqrt(sum((A*dataunif.x - b*ones(1,length(itersunif))).^2));

itersrow = datarow.iter;
errrow = sqrt(sum((datarow.x - xhat*ones(1,length(itersrow))).^2));
resrow = sqrt(sum((A*datarow.x - b*ones(1,length(itersrow))).^2));

%%
clf
semilogy(itersoptn,erroptn,itersoptl,erroptl,itersunif,errunif,itersrow,errrow)
legend({'opt $\|I-V^TDA\|$','opt $\lambda$','unif','$\propto \langle a_i,v_i\rangle$'},'Interpreter','latex')
xlabel('$k$','Interpreter','latex')
ylabel('$\|x_k-\hat x\|$','Interpreter','latex')
if writeout
    matlab2tikz('width','\figurewidth',...
        'extraaxisoptions',['legend style={font=\scriptsize},'], ...
        '../tex/figures/example_optimize_prob_error.tex');
end

%%
clf
plot(poptn),hold on, plot(poptl),plot(punif),plot(prow),hold off
legend({'opt $\|I-V^TDA\|$','opt $\lambda$','unif','$\propto \langle a_i,v_i\rangle$'},'Interpreter','latex')
xlabel('$i$','Interpreter','latex')
ylabel('$p_i$','Interpreter','latex')
if writeout
    matlab2tikz('width','\figurewidth',...
        'extraaxisoptions',['legend style={font=\scriptsize},'], ...
        '../tex/figures/example_optimize_probs.tex');
end


%% convergence rates and expected improvements
I = eye(n);
D = diag(punif./scpav);
rhounif = max(abs(eig(I - V'*D*A)));
lambdaunif = min(eig(V'*D*A + A'*D*V - A'*S*D*A));
normunif = norm(I-V'*D*A);

D = diag(prow./scpav);
rhorow = max(abs(eig(I - V'*D*A)));
lambdarow = min(eig(V'*D*A + A'*D*V - A'*S*D*A));
normrow = norm(I-V'*D*A);

D = diag(poptl./scpav);
rhooptl = max(abs(eig(I - V'*D*A)));
lambdaoptl = min(eig(V'*D*A + A'*D*V - A'*S*D*A));
normoptl = norm(I-V'*D*A);

D = diag(poptn./scpav);
rhooptn = max(abs(eig(I - V'*D*A)));
lambdaoptn = min(eig(V'*D*A + A'*D*V - A'*S*D*A));
normoptn = norm(I-V'*D*A);

fprintf('\n')
fprintf('           unif      row       optl      optn\n')
fprintf('1-lambda: %f, %f, %f, %f\n', 1-lambdaunif,1-lambdarow,1-lambdaoptl,1-lambdaoptn)
fprintf('rho:      %f, %f, %f, %f\n', rhounif,rhorow,rhooptl,rhooptn) 
fprintf('norm:     %f, %f, %f, %f\n', normunif,normrow,normoptl,normoptn) 

