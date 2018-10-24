function [Zexp,H] = binmatfac_set(X, Z, Zsets,fig_nr)

% Requires:     select_column_set.m


nZets = length(Zsets)
FPN = zeros(1, nZets);
[n,d] = size(X);
for set = 1: nZets
  Zet = Z(:,Zsets{set});
  K = size(Zet,2);
  
  Zexp = Zet;
  % Expand the set to all cominations
  for k = 2 : K-1
    c = combnk(1:K,k);
    for i = 1: size(c,1)
      Zexp = [Zexp logical(sum(Zet(:,c(i,:)),2))];
    end
  end
  Zexp = [Zexp logical(sum(Zet(:,1:K),2))];
      
  % figure(40), imagesc(Zexp), colormap(gray), 
  
  H = calculate_h(X,Zexp);
  
  A = logical(Zexp*H);
  FP = sum(sum(A & ~X));
  FN = sum(sum(~A & X));
  figure(fig_nr), imagesc(A), colormap(gray), title([FP FN])
  xlabel(set)
  FPN(set) = FP+FN;
    
%   miss = n*d;
%   W = false(n,K); 
%   H = false(K,d);
%   mask = false(size(X));
%   tic
%   for  it = 1: 10
%     Zet = Z(:,Zsets{set});
%     for k = 1: K
%       [w, h, Zet] = select_column_set(X,Zet,mask);
%       W(:,k) = w;
%       H(k,:) = h;
%       mask = mask | w*h;
%     end  
% 
%     FP = sum(sum(mask & ~X));
%     FN = sum(sum(~mask & X));
%    
%     if miss <= FP+FN
%       figure(fig_nr), imagesc(mask), colormap(gray), title([FP FN])
%       xlabel(set),  drawnow
%       break
%     else
%       miss = FP+FN;
%       FPN(set) = miss;
%     end
%   end
  toc
end  

[~,idx] = min(FPN)

 set = idx
  Zet = Z(:,Zsets{set});
  K = size(Zet,2)
  Zsets{set}
  Zexp = Zet;
  % Expand the set to all cominations
  for k = 2 : K-1
    c = combnk(1:K,k);
    for i = 1: size(c,1)
      Zexp = [Zexp logical(sum(Zet(:,c(i,:)),2))];
    end
  end
  Zexp = [Zexp logical(sum(Zet(:,1:K),2))];
      
  figure(40), imagesc(Zet), colormap(gray), 
  
  H = calculate_h(X,Zexp);
  
  A = logical(Zexp*H);
  FP = sum(sum(A & ~X));
  FN = sum(sum(~A & X));
  figure(fig_nr), imagesc(A), colormap(gray), title([FP FN])
  xlabel(set)
  FPN = FP+FN
  pause

% for set = idx
%   Zet = Z(:,Zsets{set});
%   K = size(Zet,2)
%   miss = n*d;
%   W = false(n,K); 
%   H = false(K,d);
%   W0 = false(n,K); 
%   H0 = false(K,d);
%   mask = false(size(X));
%   
%   for  it = 1: 10
%     Zet = Z(:,Zsets{set});
%     for k = 1: K
%       [w, h, Zet] = select_column_set(X,Zet,mask);
%       W(:,k) = w;
%       H(k,:) = h;
%       mask = mask | w*h;
%     end  
% 
%     FP = sum(sum(mask & ~X));
%     FN = sum(sum(~mask & X));
%     figure(fig_nr), imagesc(mask), colormap(gray), title([FP FN])
%     xlabel(idx),  drawnow
%     if miss <= FP+FN
%       break
%     else
%       miss = FP+FN;
%       W0 = W;
%       H0 = H;
%       FPN(set) = miss;
%     end
%   end
% end  
