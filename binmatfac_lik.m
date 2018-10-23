function [W,H] = binmatfac_lik(X, Z, K, min_class, max_class, fig_nr)

% Requires:     select_column_all.m

% Before the probability parameters are estimated
[n,d] = size(X);
mask = false(n,d);

% Start up: find the columns that embraces all cells, and select the best
sumZ = sum(Z,1);
midx = find(sumZ > max_class);  % These are candidates 
m = size(Z,2);                  % Number of columns in Z
H = false(m,d);                 % The row vectors that will be made

Lw = zeros(m,d);
for col = 1: m
  w = Z(:,col);
  for j = 1: d
    WndX = w == X(:,j); 
    Lw(col,j) = sum(WndX);
  end    
end
[~, idx] = max(Lw);          
for col = 1: m
  H(col,:) = idx == col;                 
end

crit = zeros(1,m);
for col = 1: m
 A = Z(:,col)*H(col,:);
 eq = A == X;
 crit(col) = sum(eq(:));
end
[~,best_col] = max(crit);     % Which w has most 1's in the h
midx(midx == best_col) = [];
Z(:,midx) = []; % Delete the others

%%%%%%%%%%%%%% Now for the real column selection %%%%%%%%%%%%%%%%%%%%%%%%%
%% 
W = false(n,K); 
H = false(K,d);
for k = 1: K
  tic
  [w, h, Z] = select_column_all(X,Z,mask, min_class);
  W(:,k) = w;
 
  H(k,:) = h;
  mask = mask | w*h;
  toc
  FP = sum(sum(mask & ~X));
  FN = sum(sum(~mask & X));
  figure(fig_nr), imagesc(mask), colormap(gray), title([FP FN])
  xlabel(k)
  drawnow
  if FP > FN % this did not work
    break
  end
end  

%% Update the row vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Z = W;
mask = false(n,d);
for k = 1: K
  if size(Z,2) < 1
    break
  end
  [w, h, Z] = select_column_all(X,Z,mask, min_class);
  W(:,k) = w;
  H(k,:) = h;
  
  mask = mask | w*h;
  figure(fig_nr), imagesc(mask), colormap(gray), title(k), drawnow
end  
