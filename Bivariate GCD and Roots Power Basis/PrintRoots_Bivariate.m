function [] = PrintRoots_Bivariate(roots_x, roots_y, polys_xy)

nRoots_x = size(roots_x,1);
nRoots_y = size(roots_y,1);

nFactors_xy = size(polys_xy,1);

% for each root in terms of x
str = sprintf('f(x) =' );

str = strcat(str,'... \n');

for i = 1:1:nRoots_x
    
    root = roots_x(i,1);
    mult = roots_x(i,2);
    
    str_a = '(x ';
    if root >=0 % root is positive so (x-r)
        str_b = sprintf(' - %2.4f)^{%i}', abs(root), mult);
    else % root is negative (x+abs(r))
        str_b = sprintf(' + %2.4f)^{%i}', abs(root), mult);
    end
    temp_str = strcat(str_a,str_b);
    
    str = strcat(str,temp_str);
end

str = strcat(str,'... \n');

for i = 1:1:nRoots_y
    
    root = roots_y(i,1);
    mult = roots_y(i,2);
    
    str_a = '(y ';
    if root >=0 % root is positive so (x-r)
        str_b = sprintf(' - %2.4f)^{%i}', abs(root), mult);
    else % root is negative (x+abs(r))
        str_b = sprintf(' + %2.4f)^{%i}', abs(root), mult);
    end
    temp_str = strcat(str_a,str_b);
    
    str = strcat(str,temp_str);
end

str = strcat(str,'... \n');

for k = 1:1:nFactors_xy
    
    coef_matrix = polys_xy{k};
    
    [r,c] = size(coef_matrix);
    m1 = r - 1;
    m2 = c - 1;
    
    % Get the number of diagonals in the matrix of coefficients of f(x,y)
    num_diags = (m1+1)+(m2+1)-1;
       
    temp_str = '(';
    
    
    for tot = 0:1:num_diags
        for i = tot:-1:0
            j = tot - i;
            
            if(i > m1 || j > m2)
            else
                
                temp_str = strcat(temp_str,...
                    sprintf('%2.2f',coef_matrix(i+1,j+1)));
                
                if i ~= 0
                    temp_str = strcat(temp_str,sprintf('x^{%i}',i));
                end
                if j ~= 0
                    temp_str = strcat(temp_str,sprintf('y^{%i}',j)); 
                end
            
                temp_str = strcat(temp_str,'+');
            end
            
            
            
        end
        
    end
    
    temp_str = strcat(temp_str, ')');
    str = strcat(str,temp_str);
end



fprintf(str);
fprintf('\n')





end