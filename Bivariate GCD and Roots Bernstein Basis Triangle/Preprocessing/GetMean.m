function [lambda] = GetMean(fxy, m, n_k)
%
% % Inputs
%
% fxy : (Matrix)
%
% m : (Int)
%
% n_k : (Int)
%
% % Outputs
%
% lambda : (Float)



global SETTINGS

switch SETTINGS.MEAN_METHOD
    
    case 'Arithmetic Mean'
        lambda = GetArithmeticMean(fxy,m,n_k);
    
    case 'Geometric Mean My Method'
        
        lambda = GeometricMean_MyMethod(fxy,m,n_k);
        
    case 'Geometric Mean Matlab Method'
        
        lambda = GeometricMean_MatlabMethod(fxy,m,n_k);
        
    case 'None'
        
        lambda = 1;
        
    otherwise
        error([mfilename ' : Geometric Mean is either (My Method) or (Matlab Method)'])
end

end
