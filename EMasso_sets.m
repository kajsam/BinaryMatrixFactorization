function [A, Z, Zsets] = EMasso_sets(X, Z0, Zsets, A0, fig_nr)

% Requires:     association_matrix_otsu.m, binmatfac_lik.m, 
%               select_column_likelihood.m 

% Input:        X - binary matrix of gene expression
%               K - maximum rank
%               A1 - Approximated structure matrix           

[n,d] = size(X);
n1 = sum(sum(X));

A = cell(1,5);
Z = cell(1,5);

min_class = ceil(0.001*n) % Minimum class size
max_class = floor(0.95*n) % Maximum class size

%% First time is always special

if Z0 == 0
    tic
  Z{1} = association_matrix_otsu(X, min_class, max_class, fig_nr);
  toc
  % figure(fig_nr), imagesc(Z{1})
  title('Candidate columns'), colormap(gray)
  % Z{1} = reduxZoverlap(Z{1}, min_class);
else
  Z{1} = Z0;      
end
drawnow

figure(fig_nr+1), imagesc(Z{1}), colormap(gray), title(size(Z{1},2))
drawnow

sumZ = sum(Z{1},1);
[~, idx] = sort(sumZ,'descend');
Z{1} = Z{1}(:,idx);
figure, imagesc(Z{1}), colormap(gray), title('Candidate columns sorted')

if isempty(Zsets)
    tic
  Zsets = column_set(Z{1},0);
  toc
  
  
end

qual = zeros(1, length(Zsets));
  for set = 1: length(Zsets)
    q = sum(Z{1}(:,Zsets{set}),2);
    qual(set) = sum(~q)+ sum(q>1);
  end
  [~, qidx] = sort(qual,'ascend');
Zsets_sort = Zsets(qidx);




if A0 == 0
  % Association matrix
  Zsets_half = Zsets_sort(1:floor(length(Zsets_sort)/2));
  tStart = tic;
  [W,H] = binmatfac_set(X, Z{1}, Zsets_half, 0);
  tElapsed = toc(tStart)
  A{1} = logical(W*H);
else
  A{1} = A0; 
end
figure(fig_nr), imagesc(Z{1}), 
title('Candidate columns'), colormap(gray)
figure(fig_nr+1), imagesc(A{1}), colormap(gray)  
title(strcat('Initial A, K = ',' ',num2str(size(W,2))))
[n1 sum(sum(A{1}))]
% Complete rounds
return
for r = 1: 2
  fig_nr = fig_nr+3;
  [alphabeta_c, cell_effect] = alpha_beta_cell(X, A{r}, 0);

  Z{r+1} = association_matrix_cell_adjusted(X, A{r}, alphabeta_c(end,:), min_class, fig_nr-1);
  % pause
  % Z{r+1} = cell_association_matrix(X, A{r}, cell_effect,[fig_nr fig_nr+1]);
  ylabel(r+1)
  %cc = size(Z{r+1},2);
  K = min(K,size(Z{r+1},2));
  if K < 3
    break
  end

  [W,H] = binmatfac_cell_lik(X, Z{r+1}, K,alphabeta_c, fig_nr+2);
  A{r+1} = logical(W*H);
  [n1 sum(sum(A{r+1}))]
  title(strcat(num2str(r),'th round A, K = ',' ',num2str(min(K,size(Z{r+1},2)))))
end

