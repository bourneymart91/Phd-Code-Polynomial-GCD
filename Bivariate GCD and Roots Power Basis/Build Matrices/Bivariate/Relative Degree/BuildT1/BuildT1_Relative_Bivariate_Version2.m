function T1 = BuildT1_Relative_Bivariate_Version2(fxy, n1_k1, n2_k2)
%
% % Inputs
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% n1_k1 : (Int) Degree of v(x,y) with respect to x
%
% n2_k2 : (Int) Degree of v(x,y) with respect to y
%
% % Outputs
%
% T1 : (Matrix) 

% Get the degree of f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get array of polynomials f_{i}(x) so that f = f0*1 + f1*y + f2*y^2 + ...
arr_fx = cell(m2+1,1);
for i = 0:1:m2
   arr_fx{i+1} = fxy(:,i+1); 
end


% Get array of partitions of the Sylvester matrix S_{k1,k2}(f(x,y),g(x,y))
arr_T1 = cell(m2 + n2_k2 + 1, n2_k2 + 1);
[arr_T1{:}] = deal(zeros(m1+n1_k1 + 1, n1_k1+1));

for j = 0:1:n2_k2
    for i = j : 1 : j + m2
        
        arr_T1{i+1,j+1} = BuildT1_Univariate(arr_fx{i-j+1}, n1_k1);
        
    end
end

T1 = cell2mat(arr_T1);

end

