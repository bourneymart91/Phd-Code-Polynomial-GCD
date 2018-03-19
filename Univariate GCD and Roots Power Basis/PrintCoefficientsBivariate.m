function [] = PrintCoefficientsBivariate(fxy,f)
% Given the polynomial f(x,y) print out the polynomial.
%
%   Input.
%
%   fxy : Matrix of coefficients of polynomial f(x,y), where (i-1,j-1)
%   entry is the coefficient of a_{i,j}x^{i}y^{j}.
%
%   f   : String indicating the name of the function 'f' or 'g'

% Get the degree of polynomial f(x,y)
[m1,m2] = GetDegree_Bivariate(fxy);

if m2 == 0
    str = sprintf('%s(x) = ',f);
else
% Initialise output string
str = sprintf('%s(x,y) = ',f);
end


% for each row of matrix f(x,y)
for i = 0:1:m1
    % for each column of matrix f(x,y)
    for j = 0:1:m2
        
        
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
                str_var_x = sprintf('x^{%i}',i);
            else
                str_var_x = '';
            end
            % Get String for powers of y
            if (j~=0)
                str_var_y = sprintf('y^{%i}',j);
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
