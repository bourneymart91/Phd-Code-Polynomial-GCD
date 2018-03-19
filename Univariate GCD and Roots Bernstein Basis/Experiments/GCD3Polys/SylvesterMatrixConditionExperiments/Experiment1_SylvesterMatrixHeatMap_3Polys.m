function [] = Experiment1_SylvesterMatrixHeatMap_3Polys(m, n, o, k)
% This experiment plots heat maps of the coefficient multipliers in the
% subresultant matrices
%
% Consider three alternative orderings of polynomials f(x), g(x) and h(x)
% in the k-th sylvester subresutlant matrix
% 
% % Inputs
%
% m : (Int) The degree of the polynomial f(x)
% 
% n : (Int) The degree of the polynomial g(x)
%
% o : (Int) The degree of the polynomial h(x) 
%
% k : (Int) Index of k-th subresultant matrix
%
%
% % Examples
%
% Experiment1_SylvesterMatrixHeatMap_3Polys(28, 18, 19, 1)

close all;
clc;

% Get heat map of \hat{S}(f(x), g(x), h(x))
SylvesterMatrixHeatMap_3Polys(m, n, o, k);

% Get heat map of \hat{S}(g(x), f(x), h(x))
SylvesterMatrixHeatMap_3Polys(n, m, o, k);

% Get heat map of S(h,f,g)
SylvesterMatrixHeatMap_3Polys(o, m, n, k);

end