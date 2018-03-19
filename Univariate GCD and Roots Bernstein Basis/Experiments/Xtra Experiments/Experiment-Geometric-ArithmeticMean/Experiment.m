function [] = Experiment(ex_num)
% An experiment where the coefficients of polynomials f(x) and g(x) are 
% plotted alongside their forms normalised by geometric and arithmetic mean
% 
% % Input
% 
% ex_num : (String) Example number
%
%

close all;
clc;

% Get the coefficients of the polynomials f(x) and g(x)
[fx, gx] = Examples_GCD_FromCoefficients(ex_num);

% Get degree of polynomials f(x) and g(x)
m = GetDegree(fx);
n = GetDegree(gx);

% Set k = 2
k = 2;

global SETTINGS
SETTINGS.BOOL_LOG = false;

% SYLVESTER MATRIX_VARIANT
%   'T'
%   'DT'
%   'TQ'
%   'DTQ'
SETTINGS.SYLVESTER_MATRIX_VARIANT = 'DTQ';

% Get Geometric mean of the entries of f(x) and g(x) in the k-th
% subresultant matrix
gm_fx = GetGeometricMeanMatlabMethod(fx, n - k);
gm_gx = GetGeometricMeanMatlabMethod(gx, m - k);


% Get Arithmetic Mean of the entries of f(x) and g(x) in the k-th
% subresultant matrix
am_fx = GetArithmeticMean(fx, n - k);
am_gx = GetArithmeticMean(gx, m - k);

% Get f(x) and g(x) normalised by their geometric means
fx_n_gm = fx ./ gm_fx;
gx_n_gm = gx ./ gm_gx;

% Get f(x) and g(x) normalised by their arithmetic means
fx_n_am = fx ./ am_fx;
gx_n_am = gx ./ am_gx;

% plot on log scale
bool_plot_logs = true;

if bool_plot_logs
    fx_n_gm = log10( fx_n_gm);
    fx_n_am = log10( fx_n_am);
    fx = log10(fx);
    
    gx_n_gm = log10( gx_n_gm);
    gx_n_am = log10( gx_n_am);
    gx = log10(gx);
end

% Plot f(x) and normalised forms of f(x)
figure_name = 'fx and fx_n';
figure('name',figure_name)
hold on
x_vec = 0:1:m;
plot(x_vec, fx_n_gm, '-s', 'DisplayName', '\bar{f}(x)_{gm}')
plot(x_vec, fx_n_am, '-s', 'DisplayName', '\bar{f}(x)_{am}')
plot(x_vec, fx, '-s', 'DisplayName', 'f(x)')
legend(gca, 'show');
hold off

% Plot g(x) and normalised forms of g(x)
figure_name = 'gx and gx_n';
figure('name', figure_name)
hold on
x_vec = 0:1:n;
plot(x_vec, gx_n_gm, '-s', 'DisplayName', '\bar{f}(x)_{gm}')
plot(x_vec, gx_n_am, '-s', 'DisplayName', '\bar{f}(x)_{am}')
plot(x_vec, gx, '-s', 'DisplayName', 'f(x)')
legend(gca, 'show');
hold off



end



