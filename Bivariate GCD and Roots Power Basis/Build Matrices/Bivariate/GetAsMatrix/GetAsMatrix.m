function fxy_matrix = GetAsMatrix (fxy_vec, m1, m2)
%
% % Inputs
%
% fxy_vec : (Vector) Coefficients of polynomial f(x,y) as a vector
%
% m1 : (Int) Degree of f(x,y) with respect to x
%
% m2 : (Int) Degree of f(x,y) with respect to y

global SETTINGS

switch SETTINGS.VECTORISATION_METHOD
    
    case 'Version 1'
        
        fxy_matrix = GetAsMatrix_Version1(fxy_vec, m1, m2);
        
    case 'Version 2'
        
        fxy_matrix = GetAsMatrix_Version2(fxy_vec, m1, m2);
        
    otherwise
        
        error('Error')
        
end
end 
