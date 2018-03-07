addpath('AIRtools1.3/')
rng(8219)
% Assemble a projection matrix
N = 50; %pixels per row and column
theta = 0:5:179; % angles
p = 150;
d = 1.4*N;
[Afull,b,x,theta,p,d] = paralleltomo(N,theta,p,d);
[mfull,n] = size(Afull);

% A by taking every third row (=ray) only
subsamp = 3;
A = Afull(ceil(subsamp/2):subsamp:end,:);
rownormsA = sqrt(sum(A.^2,2));
A = A(rownormsA~=0,:);

[m,n] = size(A);
% V by averaging 3 consecutive rows (equiv to "wider" rays)
w = ones(1,subsamp);
P = kron(speye(mfull/subsamp),w/sum(w));
V = P*Afull;
V = V(rownormsA~=0,:);

%%

%%
% orthogonal basis for the range ov V^T 
%ZV = orth(full(V'));
%ZA = orth(full(A'));

% a smooth image from AIRtools
im = phantomgallery('ppower',N,0.3,1.3,155432);
im = imfilter(im,fspecial('gaussian',16,4));
im = im/max(max(im));
x = im(:);
%find a nice x in the range of V'
% project it onto range of V^T
%x = ZV*(ZV'*x);
% actually already more or less in the range of V^T

%
subplot(1,2,1)
imagesc(reshape(x,N,N))
%subplot(1,2,2)
%imagesc(reshape(ZA*(ZA'*x),N,N))

%% Reconstruct by RKMA and RK

writeout = true;
b = A*x;
maxiter = 20*m;
av = diag(A*V');
probsV = av./sum(av);ones(m,1)/m;
norma2 = sum(A.^2,2);
probsA = norma2/sum(norma2);


[xV,dataV] = rkma(A,V,b,probsV,maxiter);
norm(xV-x)
[xA,dataA] = rkma(A,A,b,probsA,maxiter);

figure(1)
subplot(3,2,1)
imagesc(reshape(dataV.x(:,2),N,N))
subplot(3,2,2)
imagesc(reshape(dataA.x(:,2),N,N))
subplot(3,2,3)
imagesc(reshape(dataV.x(:,4),N,N))
subplot(3,2,4)
imagesc(reshape(dataA.x(:,4),N,N))
subplot(3,2,5)
imagesc(reshape(dataV.x(:,end),N,N))
subplot(3,2,6)
imagesc(reshape(dataA.x(:,end),N,N))
norm(xA-x)

if writeout
    imwrite(reshape(x,N,N),'../tex/images/ct_xhat.png')
    imwrite(reshape(dataV.x(:,2),N,N),'../tex/images/ct_rkma_1_sweep.png')
    imwrite(reshape(dataV.x(:,4),N,N),'../tex/images/ct_rkma_3_sweeps.png')
    imwrite(reshape(dataV.x(:,11),N,N),'../tex/images/ct_rkma_10_sweeps.png')
    imwrite(reshape(dataV.x(:,21),N,N),'../tex/images/ct_rkma_20_sweeps.png')
    imwrite(reshape(dataA.x(:,2),N,N),'../tex/images/ct_rk_1_sweep.png')
    imwrite(reshape(dataA.x(:,4),N,N),'../tex/images/ct_rk_3_sweeps.png')
    imwrite(reshape(dataA.x(:,11),N,N),'../tex/images/ct_rk_10_sweeps.png')
    imwrite(reshape(dataA.x(:,21),N,N),'../tex/images/ct_rk_20_sweeps.png')
end
%%
itersV = dataV.iter;
itersA = dataA.iter;
errV = sqrt(sum((dataV.x - x*ones(1,length(itersV))).^2));
errA = sqrt(sum((dataA.x - x*ones(1,length(itersA))).^2));
figure(2)
semilogy(itersV,errA,itersA,errV)
xlabel('$k$','Interpreter','latex')
ylabel('$\|x_k-\hat x\|$','Interpreter','latex')
legend('RK','RKMA')

if writeout
    matlab2tikz('width','\figurewidth',...
        'extraaxisoptions',['legend style={font=\scriptsize},'], ...
        '../tex/figures/example_ct_error.tex');
end


