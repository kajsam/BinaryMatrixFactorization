function [W,H] = binmatfac_cell_lik(X, Z, K, cell_effect)

% Requires:     select_column_cell_likelihood.m
[n,d] = size(X);
mask = false(n,d);
% mask = logical(randi([0 1],n,d));

W = false(n,K); 
H = false(K,d);
for k = 1: K
  tic
  [w, h, Z] = select_column_cell_likelihood(X,Z,mask, cell_effect);
  W(:,k) = w;
 
  H(k,:) = h;
  mask = mask | w*h;
  toc
  k
end  


