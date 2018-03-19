function T1 = BuildT1_Relative_Bivariate(fxy, n1_k1, n2_k2)
%
% % Inputs
%
% fxy : (Matrix)
%
% n1_k1 : (Int)
%
% n2_k2 : (Int)

global SETTINGS

switch SETTINGS.VECTORISATION_METHOD
    
    case 'Version 1'
        
        T1 = BuildT1_Relative_Bivariate_Version1(fxy, n1_k1, n2_k2);
        
    case 'Version 2'
        
        T1 = BuildT1_Relative_Bivariate_Version2(fxy, n1_k1, n2_k2);
        
    otherwise
        error('Error in Vectorisation Method')
end

        

end



