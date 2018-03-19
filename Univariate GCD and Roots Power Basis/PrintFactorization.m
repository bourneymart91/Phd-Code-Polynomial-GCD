function [] = PrintFactorization(fx_roots,f)

% Get number of distinct roots in f
[r,~] = size(fx_roots);
nDistinctRoots = r;

% Begin the string
str = strcat(f, '(x) = ');

% For each distinct root, add the factor (x-r)^{m} to the string
for i = 1 : 1 : nDistinctRoots
    
    % Get the root
    root = fx_roots(i,1);
    
    % Get the multiplicity
    multiplicity = fx_roots(i,2);
    
    % Append to the string
    if root < 0
        str_temp = sprintf('(x + %3.1f)^{%i}',abs(root),multiplicity);
    else % root > 0
        str_temp = sprintf('(x - %3.1f)^{%i}',root,multiplicity);
    end
    str = strcat(str,str_temp);
    
end

% New line character
str = strcat(str,'\n');

fprintf(str)

end