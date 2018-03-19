function [lambda] = GetMean_2Partitions(fxy, m, n, o, k)
%
%
% % Inputs
%
% fxy : (Matrix) Coefficients of f(x,y)
%
%
%
%


global SETTINGS

switch SETTINGS.MEAN_METHOD
    
    case 'Arithmetic Mean'
        error('Code not complete')
        lambda = GetArithmeticMean(fxy, m, n_k);
    
    case 'Geometric Mean My Method'
        error('Code not complete')
        lambda = GeometricMean_MyMethod(fxy,m,n_k);
        
    case 'Geometric Mean Matlab Method'
        
        lambda = GeometricMean_MatlabMethod_2Partitions(fxy, m, n, o, k);
        
    case 'None'
        
        lambda = 1;
        
    otherwise
        error([mfilename ' : Geometric Mean is either (My Method) or (Matlab Method)'])
end

end
