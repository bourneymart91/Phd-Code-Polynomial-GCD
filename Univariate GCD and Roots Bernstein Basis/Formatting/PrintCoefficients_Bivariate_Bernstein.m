function [] = PrintCoefficients_Bivariate_Bernstein(fxy, poly_name)
% Given the polynomial f(x,y) print out the polynomial.
%
% % Input.
%
%
% fxy : (Matrix) Matrix of the coefficients of the polynomial f(x,y), where
% (i-1,j-1) entry is the coefficient of a_{i,j}B_{i}(x)B_{j}(y).
%
% f : (String) The name of the function 'f' or 'g'



% Get the degree of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);



if m2 == 0
    str = sprintf('%s(x) = ',poly_name);
else
    
    % Initialise output string
    str = sprintf('%s(x,y) = ',poly_name);
end


% For each row of matrix of coefficients of f(x,y)
for i = 0 : 1 : m1
    
    % For each column of matrix of coefficients of f(x,y)
    for j = 0 : 1 : m2
        
        
        if (i==0 && j == 0) % If first coefficient
            str_sign = '';
            str_variables = '';
        else
            if fxy(i+1,j+1) >= 0 % If coefficient is > 0
                str_sign = ' +';
            else
                str_sign = ' ';
            end
            
            % Get String for powers of x
            if (i~=0)
                str_var_x = sprintf('B_{%i}^{%i}(x)',i,m1);
            else
                str_var_x = '';
            end
            % Get String for powers of y
            if (j~=0)
                str_var_y = sprintf('B_{%i}^{%i}(y)',j,m2);
            else
                str_var_y = '';
            end
            
            % Concatenate variable strings
            str_variables = strcat(str_var_x,str_var_y);
        end
        
        
        str_coef = sprintf('%2.4f',fxy(i+1,j+1));
        
        % Concatenate string
        temp_str = strcat(str_sign,str_coef,str_variables);
        
        % Concatenate string
        str = strcat(str,temp_str);
    end
end

% Print string to screen
fprintf(str)
fprintf('\n')




end
