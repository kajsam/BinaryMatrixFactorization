function [A, Z, alphabeta_c, Pi] = EMasso_main(X, K, A0, Z0, fig_nr)

% Requires:     passociation_matrix_otsu.m, binmatfac_lik.m, 
%               select_column_likelihood.m 

% Input:        X - binary matrix of gene expression
%               K - maximum rank
%               A1 - Approximated structure matrix           

[n,d] = size(X);


%% First time is always special

if A0 == 0
  % Association matrix
  'first go'
  if Z0 == 0
    tic
    [Z0, tau] = passociation_matrix_otsu(X, X, 0, 0);
    figure(fig_nr), imagesc(Z0), 
    title(strcat('Median tau :', num2str(median(tau)))), colormap(gray)
    toc
  end
  tic
  [W,H] = binmatfac_lik(X, Z0, min(K,size(Z0,2)));
  for k = 1: K
    Ask = logical(W(:,k)*H(k,:));
    figure(fig_nr+1), subplot(1,K,k)
    imagesc(Ask), colormap(gray), xlabel(k), drawnow
    if k == 1
      ylabel('Components')
      title('Initial A')
    end
  end
  A0 = logical(W*H);
  toc

  figure(fig_nr+2), imagesc(A0), colormap(gray)  
  title(strcat('Initial A, K = ',' ',num2str(K)))
end

A = cell(1,5);
Z = cell(1,5);

A{1} = A0;
Z{1} = Z0;
[alphabeta_c, cell_effect, Pi] = alpha_beta_cell(X, A{1}, fig_nr+3);

[Z{2}, thresh] = cell_association_matrix(X, A{1}, cell_effect,fig_nr+4);

cZc = [Z{2} ~Z{2}];

whos cZc
[W,H] = binmatfac_cell_lik(X, cZc, min(K,size(Z{2},2)),alphabeta_c);
for k = 1: K
  Ask = logical(W(:,k)*H(k,:));
  figure(fig_nr+5), subplot(1,K,k)
  imagesc(Ask), colormap(gray), xlabel(k), drawnow
  if k == 1
    ylabel('Components')
    title('First cell A')
  end
end
A{2} = logical(W*H);

figure(fig_nr+6), imagesc(A{2}), colormap(gray)  
title(strcat('First cell A, K = ',' ',num2str(K)))
