function C1 = cost_matrix_PWR(t, y, p, Lmin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matJ = cost_matrix_PWR(t, y, p, Lmni)
% matrice_cout calcule la matrice cout de fisher pour
% la segmentation du signal
%   C1(a,b) = sum_{t=a}^{t=b}[log(sigma2)+(xt-mu)^2/sigma2] 
% avec mu = beta'*r_i : un polynome d'ordre p; 
% ici beta se calcule pour chaque couple (a,b).
%
% ENTREES: 
%
%        y : signal de dim(nx1) (pour linstant cette fonction
%            n'est utilisable que pour des signaux monovaries) 
%        t : domaine temporel.
%        p : ordre de regression
%        Lmin : nbre de points minimum dans un segment (par defaut Lmin = 1)
%
% SORTIES:
%
%        C1 : matrice cout (partition en un seul segment) de dim [nxn]
%
%
%
% Faicel Chamroukhi, 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[n, d] = size(y);
if nargin<4
   Lmin = 1;
end
n = length(y);

X = designmatrix(t,p);

nl = n-Lmin+1;

C1 = Inf(n,n);
C1 = tril(C1,Lmin-2);

for a = 0:nl
    for b = a+1+Lmin:n     % ici et dans ce qui suit a+1 car en atlab les indices commencent de 1
        yab = y(a+1:b,:);
        X_ab = X(a+1:b,:);
        nk = b-a;%length(xab)
        beta = inv(X_ab'*X_ab)*X_ab'*yab;
        z = yab - X_ab*beta;
        sigma2 = z'*z/nk;
        C1(a+1,b)= nk * 0.5*log(2*pi) + nk*0.5*log(sigma2+eps) + (z'*z)/sigma2;
    end
    
end
