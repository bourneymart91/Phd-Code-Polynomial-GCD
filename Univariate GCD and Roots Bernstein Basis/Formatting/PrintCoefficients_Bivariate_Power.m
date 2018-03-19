function [] = PrintCoefficients_Bivariate_Power(fxy,f)
% Given the polynomial f(x,y) print out the polynomial.

% Get degree of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Initialise string
str = sprintf('%s(x,y) = ',f);

% for each row of coefficients f(x,y)
for i = 0:1:m1
    % for each column of coefficients f(x,y)
    for j = 0:1:m2
        
        
        if (i==0 && j == 0) % If the first coefficient
            str_sign = '';
            str_monomials = '';
           
        else % otherwise
            
            if fxy(i+1,j+1) >= 0
               str_sign = ' +';
            else 
                str_sign = ' ';
            end
            str_monomials = sprintf('x^{%i}y^{%i} ',i,j);  
        end
        
        str_coef = sprintf('%2.4f',fxy(i+1,j+1));
        
        temp_str = strcat(str_sign,str_coef,str_monomials);
                
        
        str = strcat(str,temp_str);
    end
end

fprintf(str)
fprintf('\n')
