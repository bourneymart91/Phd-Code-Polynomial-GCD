function plotMaxMinRowSum(vMaxRowNormR1, vMinRowNormR1, limits_k, limits_t, rank_range)
% Plots the max : min ratio of the row sums of R_{1, k}  for a range of k
% values, where k is the index of the k-th subresultant matrix and R_{1,k}
% is obtained by QR decomposition of S_{k}(f,g)
%
% % Inputs
%
% vMaxRowNormR1 : (Vector) Maximum of the sum of the rows of R_{1,k} for
% various values of k
%
% vMinRowNormR1 : (Vector) Minimum of the sum of the rows of R_{1,k} for
% various values of k
%
% limits_k : [Int Int] Range of k values for which the Condition numbers of
% S_{k} have been computed
%
% limits_t : [Int Int] Range of k which are candidates for the degree of
% the GCD of f(x) and g(x)
%
% rank_range : [(Float) (Float)]
% 
% % Outputs



global SETTINGS

%
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

%
x_vec = lowerLimit_k : 1 : upperLimit_k;

% Plot Graph of ratio of max : min row sum in R1 from the QR decompositions.
figure_name = sprintf('%s : Max:min Row Norms',mfilename);
figure('name',figure_name)

vRatio_MaxMin_RowNorm_R = vMinRowNormR1 ./ vMaxRowNormR1;
plot(x_vec, log10(vRatio_MaxMin_RowNorm_R),'red-s');
hold on

legend('Max:Min Row Norms of Rows in R1 from the QR decomposition of S_{k}');

title(sprintf('Max:Min Row Norms of Rows in R1 from the QR Decomposition of %s', ...
    SETTINGS.SYLVESTER_MATRIX_VARIANT));

xlim([1 upperLimit_k]);

hline(rank_range, {'r','r'});
vline(limits_t,{'b','b'});

grid on

hold off
end
