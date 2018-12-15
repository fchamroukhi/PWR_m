function PWR = fit_PWR_fisher(x, y, K, p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithme de fisher de regression par morceaux pour la parametrisation 
% des signaux de manoeuvre d'aiguillages en utilisant la programmation dynamique. 
% les parametres à estimer sont:
% 1. les temps de changement de différents phases du signal
%
% Une fois la partition estimée; on calcule les coeffecients de regression 
% associés a chaque segment ainsi que la variance du bruit sur chaque segment
% 2. les parametres de regression de chaque phase du signal
% 3. les variances du bruit additif sur chaque phase
%  La méthode d'estimation est le maximum de vraisemblance.
%
%
%
% Faicel Chamroukhi Decembre 2008.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(y,2)~=1, y=y';end
if size(x,2)~=1, x=x';end
Lmin=p+1;%1

n=length(y);

warning off

X = designmatrix(x,p);
%%% Initialisation : calcul de J_1;

tic
%%% matrice "coût"
C1 = cost_matrix_PPWR(x, y, p, Lmin);

%%% dynamic programming
[Ck, t_est] = dynamic_prog(C1, K);

gammak = [0 t_est(end,:)]; % change points

% estimation of the corresponding regression coefficients
mean_function = zeros(n,1);
for k=1:K
        i = gammak(k)+1;
        j = gammak(k+1);
        nk = j-i+1;
        yij = y(i:j);
        X_ij = X(i:j,:); 
        betak(:,k) = inv(X_ij'*X_ij)*X_ij'*yij;
        z = yij-X_ij*betak(:,k);
        sigma2k(k) = z'*z/nk;    %variances    
        mean_function(i:j) = X_ij*betak(:,k);
end
%
% classes estimees:
klas = zeros(n,1);
Zik = zeros(n,K);
for k = 1:K
    i = gammak(k)+1;
    j = gammak(k+1);
    klas(i:j) = k;
    Zik(i:j,k)=1;
end
PWR.param.betak = betak;
PWR.param.sigma2k = sigma2k;
PWR.param.gammak = gammak(2:end-1);%sans le 0 et le n
PWR.param.parameter_vector = [PWR.param.gammak(:);PWR.param.betak(:);PWR.param.sigma2k(:)];

PWR.stats.klas = klas;
PWR.stats.mean_function=mean_function;
PWR.stats.regressors = X*PWR.param.betak;
PWR.stats.Zik = Zik;
PWR.stats.objective = Ck(end);
PWR.stats.cputime = toc;

