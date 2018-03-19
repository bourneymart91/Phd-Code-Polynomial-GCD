function f_vec = GetAsVector(fxy_matrix)
%
%
% % Inputs
%
% fxy_matrix : (Matrix)
%
% % Outputs
%
% f_vec : (Vector)


global SETTINGS

switch SETTINGS.VECTORISATION_METHOD
    
    case 'Version 1'
        
        f_vec = GetAsVector_Version1(fxy_matrix);
        
    case 'Version 2'
        
        f_vec = GetAsVector_Version2(fxy_matrix);
        
    otherwise 
        error('Error on vectorisation method')
end



end
