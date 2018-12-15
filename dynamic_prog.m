function [J, t_est]= dynamic_prog(matJ, K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [J, t_est]= progdyn(matJ,Kmax)
% dynamic_prog : algorithme de programmation dynamique
%   matJ(i,j) = sum_{t=i}^{t=j}(xt-mu)^2 
% avec mu = beta'*r_i : un polynome d'ordre p; 
% ici beta se calcule pour chaque (i,j).
%
% ENTREES: 
%
%        matJ : matrice cout pour chaque couple (i,j) qui est J_1
%               matJ(i,j) = sum_{t=i}^{t=j}(xt-mu)^2 
%        Kmax : nbre maximal de classes.       
%
% SORTIES:
%
%        J : vecteur de dim [K x 1] qui contient la valeur du critere
%        correspondant à chaque valeur du nobmre de classe k qui varie de 1
%        jusqu'a K.
%
%
%        J = min_k (sum_k sum_{t=i}^{t=j}(xt-muk)^2)
%
%        t_est : temps de changement (partition) estimés (estimée) pour
%        chaque valeur de k. matrice de dim[Kmax,Kmax]
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n, ~]=size(matJ);

% I   = repmat(Inf, K, n);
I   = Inf(K, n);

t   = zeros(K-1, n);
I(1,:) = matJ(1,:);
if K>2
    for k=2:K-1
        for L=2:n
            [I(k,L),t(k-1,L)] = min(I(k-1,1:L-1) + matJ(2:L,L)' );
        end
    end
end
[I(K,n),t(K-1,n)] = min(I(K-1,1:n-1)+matJ(2:n,n)');
J=I(:,n);
% calculates the change point instants
t_est = diag(repmat(n,K,1));
for K=2:K
        for k=K-1:-1:1
            t_est(K,k) = t(k,t_est(K,k+1));
        end
end


