%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A polynomial piecewise regression model for the optimal segmentation of a
% time series with regime changes. It uses MLE for the estimation of the regression parameters
% with dynamic programming for the segmentation 
%
% by Faicel Chamroukhi Decembre 2008.
%
%% Please cite the following papers for this code:
%
% @article{chamroukhi_et_al_NN2009,
% 	Address = {Oxford, UK, UK},
% 	Author = {Chamroukhi, F. and Sam\'{e}, A. and Govaert, G. and Aknin, P.},
% 	Date-Added = {2014-10-22 20:08:41 +0000},
% 	Date-Modified = {2014-10-22 20:08:41 +0000},
% 	Journal = {Neural Networks},
% 	Number = {5-6},
% 	Pages = {593--602},
% 	Publisher = {Elsevier Science Ltd.},
% 	Title = {Time series modeling by a regression approach based on a latent process},
% 	Volume = {22},
% 	Year = {2009},
% 	url  = {https://chamroukhi.com/papers/Chamroukhi_Neural_Networks_2009.pdf}
% 	}
% 
% @INPROCEEDINGS{Chamroukhi-IJCNN-2009,
%   AUTHOR =       {Chamroukhi, F. and Sam\'e,  A. and Govaert, G. and Aknin, P.},
%   TITLE =        {A regression model with a hidden logistic process for feature extraction from time series},
%   BOOKTITLE =    {International Joint Conference on Neural Networks (IJCNN)},
%   YEAR =         {2009},
%   month = {June},
%   pages = {489--496},
%   Address = {Atlanta, GA},
%  url = {https://chamroukhi.com/papers/chamroukhi_ijcnn2009.pdf}
% }
% 
% @article{Chamroukhi-FDA-2018,
% 	Journal = {Wiley Interdisciplinary Reviews: Data Mining and Knowledge Discovery},
% 	Author = {Faicel Chamroukhi and Hien D. Nguyen},
% 	Note = {DOI: 10.1002/widm.1298.},
% 	Volume = {},
% 	Title = {Model-Based Clustering and Classification of Functional Data},
% 	Year = {2019},
% 	Month = {to appear},
% 	url =  {https://chamroukhi.com/papers/MBCC-FDA.pdf}
% 	}
%
%
%
% Faicel Chamroukhi Decembre 2008.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close all

%% toy time series with regime changes
% y =[randn(100,1); 7+randn(120,1);4+randn(200,1); -2+randn(100,1); 3.5+randn(150,1);]';
% n = length(y);
% x = linspace(0,1,n);

load simulated_time_series;

% model specification
K = 5;% number of segments
p = 3; % polynomial degree


pwr = fit_PWR(x, y, K, p);

fprintf('elapsed time = %g\n', pwr.stats.cputime);
%fprintf('objective value = %f\n',pwr.stats.objective);

show_PWR_results(x, y, pwr);

%% some real time series with regime changes

%load real_time_series_1
load real_time_series_2

pwr = fit_PWR(x, y, K, p);

fprintf('elapsed time = %g\n', pwr.stats.cputime);
%fprintf('objective value = %f\n',pwr.stats.objective);
yaxislim = [240, 600];
show_PWR_results(x, y, pwr, yaxislim);



