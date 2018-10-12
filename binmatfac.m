function [W,H] = binmatfac(X, Z, K)

% Requires:     select_column.m
[n,d] = size(X);
mask = false(n,d);
% mask = logical(randi([0 1],n,d));

W = false(n,K); 
H = false(K,d);
for k = 1: K
  [w, h, Z] = select_column(X,Z,mask,K-k+1);
  W(:,k) = w;
 
  H(k,:) = h;
  mask = mask | w*h;
end  


