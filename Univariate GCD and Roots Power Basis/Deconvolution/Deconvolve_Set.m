function [arr_hx] = Deconvolve_Set(arr_fx, DECONVOLUTION_METHOD)
% Deconvolve_Set
% Deconvolve the set of polynomials f_{i}(x), where the polynomails
% f_{i}(x) are outputs of a sequence of GCD computations in the Tobey and
% Horowitz algorithm.
%
% Inputs.
%
%
% arr_fx : Array of polynomials f_{i}(x)
%
% DECONVOLUTION_METHOD : (String)
%
% Outputs
%
%
% arr_hx : Array of polynomials h_{i}(x)

% Get the number of polynomials in array of f_{i}
nPolys_fx = size(arr_fx, 1);

fprintf(['Deconvolution Method : ' DECONVOLUTION_METHOD '']);

switch DECONVOLUTION_METHOD
    
    case 'Separate'
        
        % Initialise an array 
        arr_hx = cell(nPolys_fx-1,1);
        
        % For each pair f_{i}, f_{i+1} perform deconvolution
        for i = 1:1:nPolys_fx-1
            arr_hx{i} = Deconvolve(arr_fx{i}, arr_fx{i+1}) ;
        end
        
    case 'Batch'
        
        % Normalise all entries
        for i = 1:1:nPolys_fx
            arr_fx{i} = arr_fx{i}./arr_fx{i}(1,1);
        end
        
        % Deconvolve batch
        arr_hx = Deconvolve_Batch(arr_fx);
        
    case 'Batch With STLN'
        
     
        arr_hx = Deconvolve_Batch_With_STLN(arr_fx);
        
    case 'Batch Constrained'
        
        % Normalise all entries
        for i = 1:1:nPolys_fx
            arr_fx{i} = arr_fx{i}./arr_fx{i}(1,1);
        end
        
        
        % Get the degree of polynomials f_{i}(x)
        vDegt_fx = zeros(nPolys_fx,1);
        for i = 1:1:nPolys_fx
            vDegt_fx(i) = GetDegree(arr_fx{i});  
        end
                
        % Get the degree structure of the polynomials h_{i}
        vDeg_hx = diff(vDegt_fx);
        
        % Get the degree structure of the polynomials w_{i}
        vDeg_wx = diff([vDeg_hx; 0]);
        
        % Get the vector of multiplicities of each of the factors of f(x)
        vMult = find(vDeg_wx~=0);
        
        arr_hx = Deconvolve_Batch_Constrained(arr_fx, vMult);
        
    case 'Batch Constrained With STLN'
        
        % Normalise all entries
        for i = 1:1:nPolys_fx
            arr_fx{i} = arr_fx{i}./arr_fx{i}(1,1);
        end
        
        % Get the degree of polynomials f_{i}(x)
        vDegt_fx = zeros(nPolys_fx,1);
        for i = 1:1:nPolys_fx
            vDegt_fx(i) = GetDegree(arr_fx{i});  
        end
                
        % Get the degree structure of the polynomials h_{i}
        vDeg_hx = diff(vDegt_fx);
        
        % Get the degree structure of the polynomials w_{i}
        vDeg_wx = diff([vDeg_hx; 0]);
        
        % Get the vector of multiplicities of each of the factors of f(x)
        vMult = find(vDeg_wx~=0);
        
        arr_hx = Deconvolve_Batch_Constrained_With_STLN(arr_fx,vMult);
        
    otherwise
        error([...
            'SETTINGS.DECONVOLVE_METHOD'...
            ' *Separate'...
            ' *Batch'...
            ' *Batch With STLN'...
            ' *Batch Constrained'...
            ' *Batch Constrained With STLN']);
end
end


