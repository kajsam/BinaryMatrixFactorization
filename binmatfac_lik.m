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
  FP = sum(sum(mask & ~X));
  FN = sum(sum(~mask & X));
  figure(fig_nr), imagesc(mask), colormap(gray), title([FP FN]), drawnow
  xlabel(k)
  if FP > FN % this did not work
    break
  end
end  

Z = W;
mask = false(n,d);
for k = 1: K
  [w, h, Z] = select_column_all(X,Z,mask);
  W(:,k) = w;
  H(k,:) = h;
  
  mask = mask | w*h;
  figure(fig_nr), imagesc(mask), colormap(gray), title(k), drawnow
end  
