function Cf = BuildSylvesterMatrixPartition_2Polys(fxy, m, n_k)
% BuildSylvesterMatrix(fxy,gxy,m,n,k)
%
% Build the Sylvester subresultant matrix S_{k}(f,g)
% 
%
% Inputs
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% n : (Int) Total degree of polynomial g(x,y)
%
% k : (Int) Index of the Sylvester subresultant matrix S_{k} to be constructed.
%
% % Outputs
%
% Sk : (Matrix) kth Sylvester subresultant matrix S_{k}(f,g)



global SETTINGS
switch SETTINGS.SYLVESTER_MATRIX_VARIANT
    case 'T'

        % Build the matrix T_{n-k}(f(x,y))
        T1_fx = BuildT1(fxy, m, n_k);

        Cf = T1_fx;
        
    case 'DT'
        % Build the diagonal matrix D^{-1}
        D = BuildD_2Polys(m, n_k);

        % Build the matrix T_{n-k}(f(x,y))
        T1_fx = BuildT1(fxy, m, n_k);
       
        Cf = D * T1_fx;
        
    case 'DTQ'
        % Build the diagonal matrix D^{-1}
        D = BuildD_2Polys(m, n_k);

        % Build the matrix T_{n-k}(f(x,y))
        T1_fx = BuildT1(fxy, m, n_k);

        % Build Matrix Q1
        Q1 = BuildQ1(n_k);
        
        Cf = D * T1_fx * Q1;
    
    case 'TQ'

        % Build the matrix T_{n-k}(f(x,y))
        T1_fx = BuildT1(fxy, m, n_k);

        Q1 = BuildQ1(n_k);
        
        Cf = T1_fx * Q1;
        
    case 'DTQ Denominator Removed'
        
        % Build the diagonal matrix D^{-1}
        D = BuildD_2Polys(m, n_k);

        % Build the matrix T_{n-k}(f(x,y))
        T1_fx = BuildT1(fxy, m, n_k) ./ nchoosek(m + n_k, n_k);

        % Build Q1
        Q1 = BuildQ1(n_k);
        
        Cf = D * T1_fx * Q1;
        
    otherwise
        
        error('err')
        
end


end
