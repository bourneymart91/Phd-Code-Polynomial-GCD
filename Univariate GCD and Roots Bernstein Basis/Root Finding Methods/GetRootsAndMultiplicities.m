function [root_mult_matrix]= GetRootsAndMultiplicities(fx, str_method)
%
% % Inputs
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% str_method : (String) Factorisation method name
%
% % Outputs
%
% root_mult_matrix : (Matrix) where the columns contain the roots and
% corresponding multiplicity of the root.

switch str_method
    
        
    case 'Musser Method'
        root_mult_matrix = o_roots_Musser(fx);
        
    case 'Yun Method'
        root_mult_matrix = o_roots_Yun(fx);
        
    case 'Matlab Method'
        root_mult_matrix = o_roots_Matlab(fx);
        
    case 'Zeng Method'
        root_mult_matrix = o_roots_multroot(fx);
        
        
    case 'Bisection Method'
        root_mult_matrix = o_roots_Bisection(fx);
        
    case 'Subdivision Method'
        root_mult_matrix = o_roots_Subdivision(fx);
        
    case 'Bezier Clipping Method'
        root_mult_matrix = o_roots_BezierClipping(fx);
        
        
    otherwise
        error('err')
end

% Print the calculated roots and the corresponding multiplicities.
try
    PrintoutRoots(str_method, root_mult_matrix);
catch
    fprintf('No Roots Found in %s \n', str_method)
end

end