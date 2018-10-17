function [W,H] = binmatfac_lik(X, Z, K, fig_nr)

% Requires:     select_column_all.m

% Before the probability parameters are estimated
[n,d] = size(X);
mask = false(n,d);

W = false(n,K); 
H = false(K,d);
for k = 1: K
  tic
  [w, h, Z] = select_column_all(X,Z,mask);
  W(:,k) = w;
 
  H(k,:) = h;
  mask = mask | w*h;
  toc
  figure(fig_nr), imagesc(mask), colormap(gray), title(k), drawnow
end  

Z = W;
mask = false(n,d);
for k = 1: K
  tic
  [w, h, Z] = select_column_all(X,Z,mask);
  W(:,k) = w;
 
  H(k,:) = h;
  mask = mask | w*h;
  toc
  figure(fig_nr+1), imagesc(mask), colormap(gray), title(k), drawnow
end  
