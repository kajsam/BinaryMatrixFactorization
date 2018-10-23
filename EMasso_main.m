function [A, Z] = EMasso_main(X, K, Z0, A0, fig_nr)

% Requires:     passociation_matrix_otsu.m, binmatfac_lik.m, 
%               select_column_likelihood.m 

% Input:        X - binary matrix of gene expression
%               K - maximum rank
%               A1 - Approximated structure matrix           

[n,d] = size(X);
n1 = sum(sum(X));

A = cell(1,5);
Z = cell(1,5);

min_class = ceil(0.01*n) % Minimum class size
max_class = floor(0.95*n) % Maximum class size

%% First time is always special

if Z0 == 0
    tic
  Z{1} = association_matrix_otsu(X, min_class, fig_nr);
  toc
  % figure(fig_nr), imagesc(Z{1})
  title('Candidate columns'), colormap(gray)
else
  Z{1} = Z0;      
end
drawnow
Z{1} = reduxZoverlap(Z{1}, min_class);
figure(fig_nr+1), imagesc(Z{1}), colormap(gray), title(size(Z{1},2))
drawnow

K = min(K,size(Z{1},2));
if A0 == 0
  % Association matrix
  tStart = tic;
  [W,H] = binmatfac_lik(X, Z{1}, K, min_class, max_class, fig_nr+2);
  tElapsed = toc(tStart)
  A{1} = logical(W*H);
else
  A{1} = A0; 
end
figure(fig_nr), imagesc(Z{1}), 
title('Candidate columns'), colormap(gray)
figure(fig_nr+1), imagesc(A{1}), colormap(gray)  
title(strcat('Initial A, K = ',' ',num2str(K)))
[n1 sum(sum(A{1}))]
% Complete rounds
return
for r = 1: 2
  fig_nr = fig_nr+3;
  [alphabeta_c, cell_effect] = alpha_beta_cell(X, A{r}, 0);

  Z{r+1} = association_matrix_cell_adjusted(X, A{r}, alphabeta_c(end,:), min_class, fig_nr-1);
  pause
  Z{r+1} = cell_association_matrix(X, A{r}, cell_effect,[fig_nr fig_nr+1]);
  ylabel(r+1)
  cc = size(Z{r+1},2);
  K = min(K,floor(cc/2));
  if K < 3
    break
  end

  [W,H] = binmatfac_cell_lik(X, Z{r+1}, K,alphabeta_c, fig_nr+2);
  A{r+1} = logical(W*H);
  [n1 sum(sum(A{r+1}))]
  title(strcat(num2str(r),'th round A, K = ',' ',num2str(min(K,size(Z{r+1},2)))))
end

