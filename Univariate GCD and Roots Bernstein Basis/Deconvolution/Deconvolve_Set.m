function arr_hx = Deconvolve_Set(arr_fx, DECONVOLUTION_METHOD)
% Performs a series of n-1 deconvolutions over a set of n polynomials,
% where each polynomial f_{i}(x) appears in two deconvolutions, and
% h_{i}(x) = f_{i}(x)/f_{i+1}(x)
%
%
% % Inputs
%
% arr_fx : (Array of Vectors) Each cell of the array contains coefficients 
% of the polynomials f_{i}(x) 
%
% % Outputs
%
% arr_hx : (Array of Vectors) Each cell of the array contains coefficients
% of the polynomial h_{i}(x) = f_{i}(x)/f_{i+1}(x)




% Set Deconvolution Method
%   Separate : Use standard non batch method
%   Batch  : Use batch deconvolution
%   Batch With STLN
%   Batch Constrained
%   Batch Constrained With STLN
%


switch DECONVOLUTION_METHOD
    case 'Separate'
        
        % Deconvolve independent method
        arr_hx = Deconvolve_Separate(arr_fx);
        
    case 'Batch'
        
        % Deconvolve Batch Method
        arr_hx = Deconvolve_Batch(arr_fx);
        
    case 'Batch with STLN'
        
        arr_hx = Deconvolve_Batch_With_STLN(arr_fx);
        
    case 'Batch Constrained'
        
        % Get number of polynomials in batch
        nPolynomials_fx = size(arr_fx,1);
        
        % Get the degree of polynomials f_{i}(x)
        vDegree_fx = zeros(nPolynomials_fx,1);
        
        for i = 1:1:nPolynomials_fx
            vDegree_fx(i) = GetDegree(arr_fx{i});  
        end
                
        % Get the degree structure of the polynomials h_{i}
        vDegree_hx = diff(vDegree_fx);
        
        % Get the degree structure of the polynomials w_{i}
        vDegree_wx = diff([vDegree_hx; 0]);
        
        % Get the vector of multiplicity structur of factors of f_{0}(x)
        vMultiplicity = find(vDegree_wx~=0);
        
        % Get array of polynomials h_{i}(x)
        arr_hx = Deconvolve_Batch_Constrained(arr_fx, vMultiplicity);
        
    case 'Batch Constrained with STLN'
        
        % Get number of polynomials in batch
        nPolynomials_fx = size(arr_fx,1);
        
        % Get the degree of polynomials f_{i}(x)
        vDegree_fx = zeros(nPolynomials_fx,1);
        for i = 1:1:nPolynomials_fx
            vDegree_fx(i) = GetDegree(arr_fx{i});  
        end
                
        % Get the degree structure of the polynomials h_{i}
        vDegree_hx = diff(vDegree_fx);
        
        % Get the degree structure of the polynomials w_{i}
        vDegree_wx = diff([vDegree_hx; 0]);
        
        % Get vector of multiplicity of each factor in f_{0}(x)
        vMultiplicity = find(vDegree_wx~=0);
        
        % Compute the polynomials h_{i}(x)
        arr_hx = Deconvolve_Batch_Constrained_With_STLN(arr_fx,vMultiplicity);
        
    otherwise
        err_msg = sprintf(...
            [
            'SETTINGS.DECONVOLUTION_METHOD must be one of the following:\n'...
            '\t*Separate \n '...
            '\t*Batch \n '...
            '\t*Batch with STLN \n '...
            '\t*Batch Constrained\n '...
            '\t*Batch Constrained with STLN\n'...
            ]);
        error(err_msg);
end




end
