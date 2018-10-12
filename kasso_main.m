function A = kasso_main(X, K, tau, ass, A1, Z, fig_nr)

% Input:    X - binary matrix of gene expression
%           K - maximum rank
%           tau - threshold

addpath('C:\Users\kajsam\Documents\MATLAB\mdl4bmf')
addpath('C:\Users\kajsam\Documents\MATLAB\Simulation')


[n,d] = size(X)
IM = ones(n,d);
imK = max(K,3);
if ass
  'Asso'
  tic
  [W, H] = asso(double(X'), K, tau); % notice the transpose of A!
  toc
 
  W = W';
  H = H';
  As = logical(W*H);
  
  clf(figure(fig_nr))
  figure(fig_nr), subplot(3,imK,1),imagesc(IM-X), colormap(gray), title('X')
  subplot(3,imK,2), imagesc(IM-As), title(strcat('Asso cover K = ',' ',num2str(K)))
  drawnow
  
  figure(fig_nr)
  for k = 1: K
    Ask = logical(W(:,k)*H(k,:));
    subplot(3,imK,imK+k),imagesc(IM-Ask), colormap(gray), xlabel(k)
    if k == 1
      ylabel('Asso components')
    end
  end
end

%% First time is always special

if A1 == 0
  % Association matrix
  'first go'
  if Z == 0
    tic
    [Z, tau] = passociation_matrix_otsu(X, X, 0);
    figure(fig_nr+1), imagesc(Z), 
    title(strcat('Median tau :', num2str(median(tau)))), colormap(gray)
    toc
  end

  tic
  [W,H] = binmatfac_lik(X, Z, min(K,size(Z,2)));
  
  for k = 1: K
    Ask = logical(W(:,k)*H(k,:));
    figure(fig_nr), subplot(3,imK,2*imK+k)
    imagesc(IM-Ask), colormap(gray), xlabel(k), drawnow
    if k == 1
      ylabel('Uncover components')
    end
  end
  A1 = logical(W*H);
  toc

  figure(fig_nr),subplot(3,imK,3), imagesc(IM-A1) 
  title(strcat('Uncover K = ',' ',num2str(K)))
end

figure(1), subplot(1,3,2), imagesc(IM-A1), colormap(gray) 
% title(sum(sum(xor(A1,S))))
oscillations = 10;

rK = 2;
imK = ceil(oscillations/rK);
imk = 0;


A = cell(1,oscillations+1);
A{1} = A1;

fig_nr = fig_nr+2;

for osc = 1: oscillations
    osc
  % Association matrix
  tic
  Z = passociation_matrix_otsu(X, A{osc}, 0);
  Z(:,sum(Z,1)==0) =[];
  size(Z)
  toc
  imk = imk + 1;
  figure(fig_nr+1), subplot(rK, imK, imk),
  imagesc(Z), title(size(Z,2)), colormap(gray)
 
  tic
  [W,H] = binmatfac_lik(X, Z, min(K,size(Z,2)));
  A{osc+1} = logical(W*H);
  toc
  figure(fig_nr), subplot(rK, imK, imk), imagesc(IM-A{osc+1}) 
  title(osc), colormap(gray)
  drawnow
end

figure(1), subplot(1,3,3), imagesc(IM-A{osc+1}), colormap(gray)
% title(sum(sum(xor(A{osc+1},S))))

