function [fxy_noisy,noise_matrix] = AddVariableNoiseToPoly(fxy, el, eu)
%
% Add noise to the coefficients of polynomial f(x,y)
%
%
%
% Inputs
%
% fxy : matrix of coefficients of polynomial f(x,y)
%
% el : signal to noise low limit
%
% eu : signal to noise upper limit


global SETTINGS
rng(SETTINGS.SEED)

% get the degree of input polynomial f
[m1, m2] = GetDegree_Bivariate(fxy);

switch nargin
    case 2 % Only one noise is specified, set upper = lower
        
        
        rp = (2*rand(m1+1,m2+1))-ones(m1+1,m2+1);
        s = rp*el;
        
        noise_matrix = fxy.*s;
        fxy_noisy = fxy + noise_matrix;
        
        
    case 3 % Specified upper and lower bound of noise
        
        y = (2*rand(m1+1, m2+1)) - ones(m1+1, m2+1);
        s = eu *ones(m1+1, m2+1) -  y.*(eu-el);
        noise_matrix = fxy.*s;
        fxy_noisy = fxy + noise_matrix;
end