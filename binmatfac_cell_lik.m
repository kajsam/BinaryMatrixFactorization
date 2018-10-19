function [W,H] = binmatfac_cell_lik(X, Z, K, alphabeta_c, fig_nr)

% Requires:     select_column_cell_likelihood.m
[n,d] = size(X);
mask = false(n,d);
% mask = logical(randi([0 1],n,d));

W = false(n,K); 
H = false(K,d);
for k = 1: K
  [w, h, Z] = select_column_cell_likelihood(X,Z,mask, alphabeta_c);
  W(:,k) = w;
 
  H(k,:) = h;
  mask = mask | w*h;
  % figure(fig_nr), imagesc(mask), colormap(gray), title(k), drawnow
end  

Z = W;
mask = false(n,d);
for k = 1: K
  [w, h, Z] = select_column_cell_likelihood(X,Z,mask, alphabeta_c);
  W(:,k) = w;
  H(k,:) = h;
  
  mask = mask | w*h;
  figure(fig_nr), imagesc(mask), colormap(gray), title(k), drawnow
end  


