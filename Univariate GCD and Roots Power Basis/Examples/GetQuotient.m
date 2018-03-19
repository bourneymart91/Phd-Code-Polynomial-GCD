function u_roots = GetQuotient(f_roots,d_roots)
% Get the roots of quotient polynomial u(x) and corresponding multiplicities
% given the roots of polynomial f(x), and the roots of polynomial d, 
% where d is the GCD of f(x) and g(x).

% Get the number of distinct roots in f
nDistinctRoots_fx = size(f_roots,1);
u_roots = [];

% Catch the case where the exact GCD is zero, and therefore the quotient
% polynomial u is equal to f.
[r,~] = size(d_roots);
if r == 0
    u_roots = f_roots;
    return
end


% For each of the distinct roots in f(x)
for i = 1:1:nDistinctRoots_fx
    
    % Get the root
    root = f_roots(i,1);
    
    % Get corresponding multiplicity.
    mult_f = f_roots(i,2);
    
    % Look for the root in roots of d(x)
    if ~isempty(find(d_roots(:,1) == root));
        
        % Get the index of the row on which the root is found.
        [row_d,~] = find(d_roots(:,1) == root);
        
        % get the multiplicity of the root in d(x)
        mult_d = d_roots(row_d,2);
        
        % Subtract multiplicty of root in d(x) from multiplicity of the 
        % same root in f(x) to obtain multiplicity in quotient polynomial 
        % u(x)
        mult_u = mult_f - mult_d;
        
        % Add the root and its multiplicity to the set of roots for
        % quotient polynomial u(x)
        if mult_u > 0
            u_roots = [u_roots; root mult_u];
        end
        
    else
        % Otherwise, root does not exist in d(x) so exists entirely in
        % u(x).
        u_roots = [u_roots; root mult_f];
    end
end


end