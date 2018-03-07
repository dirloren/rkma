function k = sampling(p)
%function k = sampling(p)
% Return a random integer sampled with probabilities according to p, i.e.:
% p is vector of length m with non-negative entries that sum to one and the
% result of sampling(p) is k (1<=k<=m) with probability p(k).

P = cumsum(p);
k = nnz(rand>P)+1;
