function [W,H, Z, Zsets] = EMasso_cellgene(X, Z0, Zsets, W0, H0, fig_nr)

% Requires:     association_matrix_otsu.m, binmatfac_lik.m, 
%               select_column_likelihood.m 

% Input:        X - binary matrix of gene expression
%               K - maximum rank
%               A1 - Approximated structure matrix           

[n,d] = size(X);
n1 = sum(sum(X));

W = cell(1,5);
H = cell(1,5);
notWexp = cell(1,5);
Wexp = cell(1,5);
Hexp = cell(1,5);
notHexp = cell(1,5);
A = cell(1,5);
Z = cell(1,5);

min_class = ceil(0.01*n) % Minimum class size
max_class = floor(0.99*n) % Maximum class size

%% First time is always special

if Z0 == 0
    tic
  Z{1} = association_matrix_otsu(X, min_class, max_class, fig_nr);
  toc
  % figure(fig_nr), imagesc(Z{1})
  % title('Candidate columns'), colormap(gray)
  Z{1} = reduxZoverlap(Z{1}, min_class); 
  % complicating things. 
else
  Z{1} = Z0;      
end
drawnow

figure(fig_nr+1), imagesc(Z{1}), colormap(gray), title(size(Z{1},2))
drawnow

sumZ = sum(Z{1},1);
[~, idx] = sort(sumZ,'descend');
Z{1} = Z{1}(:,idx);
figure(fig_nr+2), imagesc(Z{1}), colormap(gray), title('Candidate columns sorted')
drawnow
if isempty(Zsets)
  Zsets = column_set(Z{1},0,4);
end


qual = zeros(1, length(Zsets));
for set = 1: length(Zsets)
  q = sum(Z{1}(:,Zsets{set}),2);
  qual(set) = sum(~q)+ sum(q>1);
end
[~, qidx] = sort(qual,'ascend');
Zsets_sort = Zsets(qidx);
Zsets_half = Zsets_sort(1:floor(length(Zsets_sort)/4));

if W0 == 0
  % Association matrix
  tStart = tic;
  [W{1},H{1}] = binmatfac_set(X, Z{1}, Zsets, fig_nr+2);
  
  size(W{1})
  size(H{1})
  tElapsed = toc(tStart)
  A{1} = logical(W{1}*H{1});
else 
  W{1} = W0; 
  H{1} = H0;
  A{1} = logical(W{1}*H{1});
end
figure(fig_nr), imagesc(Z{1}), 
title('Candidate columns'), colormap(gray)

[Wexp{1}, Hexp{1}, notWexp{1}, notHexp{1}] = expand(W{1},X,fig_nr);
A{1} = logical(Wexp{1}*Hexp{1});
[n1 sum(sum(A{1}))]
% Complete rounds

for r = 1: 10
  fig_nr = fig_nr +4
  r
  [pi_g] = alpha_beta_cellgene(X, Wexp{r}, Hexp{r}, notWexp{r}, notHexp{r}, 0);
   
  Z{r+1} = association_matrix_gene_adjusted(X, A{r}, pi_g(2,:), min_class, 0);
  Z{r+1} = reduxZoverlap(Z{r+1}, min_class); 
  Zsets = column_set(Z{r+1},0,0);
    
  [W{r+1},H{r+1}] = binmatfac_gene_set(X, Z{r+1}, Zsets, pi_g, fig_nr+2);
  
  [Wexp{r+1}, Hexp{r+1}, notWexp{r+1}, notHexp{r+1}] = expand_pi(W{r+1},X, pi_g, 0);
  A{r+1} = logical(Wexp{r+1}*Hexp{r+1});
  figure(fig_nr+20), subplot(1,2,1), imagesc(A{r+1}), colormap(gray)
  subplot(1,2,2), imagesc(W{r+1}), colormap(gray), title(r)

end

