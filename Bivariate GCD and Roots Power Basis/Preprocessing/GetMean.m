function lambda = GetMean(fxy)
% Get the mean of the entries of f(x,y) in the matrix C(f), by a predefined
% mean method.
%
% % Inputs
% 
% fxy : Coefficients of the polynomial f(x,y).
%
% % Outputs
%
% lamdba : Geometric mean of non-zero coefficients of polynomial f(x,y).
%
% % Note that the geometric mean is independent of the type of Sylvester
% matrix used.

global SETTINGS

switch SETTINGS.MEAN_METHOD
    case 'Geometric Mean Matlab Method'
        
        lambda  = geomean(abs(nonzeros(fxy)));
        
    case 'None'
        lambda = 1;
        
    otherwise
        error('error')
end
