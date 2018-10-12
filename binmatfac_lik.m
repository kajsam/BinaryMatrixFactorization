function [W,H] = binmatfac_lik(X, Z, K)

% Requires:     select_column.m
[n,d] = size(X);
mask = false(n,d);
% mask = logical(randi([0 1],n,d));

W = false(n,K); 
H = false(K,d);
for k = 1: K
  tic
  [w, h, Z] = select_column_likelihood(X,Z,mask);
  W(:,k) = w;
 
  H(k,:) = h;
  mask = mask | w*h;
  toc
  k
end  


