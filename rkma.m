function [x,data] = rkma(A,V,b,p,maxiter)
%function rkma(A,V,b,p)
% Perform the randomized Kaczmarz with mismatched adjoint for the equation
% Ax=b using V' instead of A' as an adjoint. The i-th row is chosen with
% probability p(i).

% pretanspose for quicker row-extraction
AT = A';
VT = V';

x = zeros(size(A,2),1);
data.iter = zeros(1,floor(maxiter/size(A,1)));
data.iter(1) = 0;
data.x = zeros(size(A,2),floor(maxiter/size(A,1)));
data.x(:,1) = x;
l = 1;
for k=1:maxiter
    i = sampling(p);
    ai = AT(:,i);
    vi = VT(:,i);
    bi = b(i);
    
    x = x - (ai'*x-bi)/(ai'*vi)*vi;
    if mod(k,size(A,1)) == 0
        l = l+1;
        data.iter(l) = k;
        data.x(:,l) = x;

        fprintf('iter = %d, residual = %e\n',k,norm(A*x-b));
    end
end
    