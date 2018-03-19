function arrRoots_d = GetDivisor(arrRoots_f,arrRoots_g)
% Given two sets of roots for polynomials f and g,
% Compare the sets of roots and return the common roots.
% 
% Inputs.
% 
% arrRoots_f : Array of roots of f(x) and corresponding multiplicities.
% 
% arrRoots_g : Array of roots of g(x) and corresponding multiplicities.
% 
% Outputs.
%
% arrRoots_d : Array of roots of d(x) and corresponding multiplicities.

% Get the number of roots in polynomial f(x,y)
num_roots_f = size(arrRoots_f,1);

% Initialise an empty array of roots for d(x,y)
arrRoots_d = [];

% for each root in f_x, check to see if it exists in g_x
for i = 1:1:num_roots_f
    
    % Get the root
    root = arrRoots_f(i,1);
    
    % Get the multiplicity of the root in f(x)
    mult_root_in_f = arrRoots_f(i,2);
    
    % Look for the root in g(x)
    [r,~] = size(arrRoots_g);
    if  r == 0
        return
    end
    
    % Find location of the root in g(x) roots array
    if ~isempty(find(arrRoots_g(:,1) == root));
        
        % Get the row index of the root in arrRoots
        [row_d,~] = find(arrRoots_g(:,1) == root);
        
        % Get the multiplicity of the root in g(x)
        mult_root_in_g = arrRoots_g(row_d,2);
        
        % Calculate the multiplicity of the root in d(x)
        mult_root_in_d = min(mult_root_in_f,mult_root_in_g);
        
        % Add the root to the array of roots for d(x)
        arrRoots_d = [arrRoots_d ; root mult_root_in_d]; 
    end
end


end

