function Sk = BuildSylvesterMatrix_3Polys(fxy, gxy, hxy, m, n, o, k)
% BuildSylvesterMatrix_3Polys(fxy, fxy2, gxy, hxy, m, n, k)
%
% Build the Sylvester subresultant matrix S_{k}(f,g)
%
%
% Inputs
%
% fxy : (Matrix)
%
% fxy2 : (Matrix)
%
% gxy : (Matrix)
%
% hxy : (Matrix) Coefficients of the polynomial f(x,y), g(x,y) and
% h(x,y)
%
% m : (Int) Degree of f(x,y)
%
% n : (Int) Degree of g(x,y)
%
% o : (Int) degree of h(x,y)
%
% k : Index of the Sylvester subresultant matrix S_{k} to be constructed.


global SETTINGS

switch SETTINGS.SYLVESTER_MATRIX_3POLY_N_EQUATIONS
    
    case '2'
        Sk = BuildSylvesterMatrix_3Polys_2Eqns(fxy, gxy, hxy, m, n, o, k);
        
    case '3'
        Sk = BuildSylvesterMatrix_3Polys_3Eqns(fxy, gxy, hxy, m, n, o, k);
        
    otherwise
        error('err');
end
end


function Sk = BuildSylvesterMatrix_3Polys_2Eqns(fxy, gxy, hxy, m, n, o, k)
% BuildSylvesterMatrix_3Polys(fxy, fxy2, gxy, hxy, m, n, k)
%
% Build the Sylvester subresultant matrix S_{k}(f,g)
%
%
% Inputs
%
% fxy : (Matrix)
%
% fxy2 : (Matrix)
%
% gxy : (Matrix)
%
% hxy : (Matrix) Coefficients of the polynomial f(x,y), g(x,y) and
% h(x,y)
%
% m : (Int) Degree of f(x,y)
%
% n : (Int) Degree of g(x,y)
%
% o : (Int) degree of h(x,y)
%
% k : Index of the Sylvester subresultant matrix S_{k} to be constructed.


global SETTINGS
switch SETTINGS.SYLVESTER_MATRIX_VARIANT
    case 'T'
        
        
        T = BuildT_3Polys_2Eqns(fxy, gxy, hxy, m, n, o, k);
        
        Sk = T;
        
    case 'DT'
        
        D = BuildD_3Polys_2Eqns(m, n, o, k);
        T = BuildT_3Polys_2Eqns(fxy, gxy, hxy, m, n, o, k);
        
        % Build the matrix DT
        Sk = D*T;
        
    case 'DTQ'
        
        
        D = BuildD_3Polys_2Eqns(m, n, o, k);
        T = BuildT_3Polys_2Eqns(fxy, gxy, hxy, m, n, o, k);
        Q = BuildQ_3Polys(m, n, o, k);
        
        % Build the matrix DTQ
        Sk = D*T*Q;
        
        
    case 'TQ'
        
        T = BuildT_3Polys_2Eqns(fxy, gxy, hxy, m, n, o, k);
        Q = BuildQ_3Polys(m, n, o, k);
        
        
        Sk = T * Q;
        
    case 'DTQ Denominator Removed'
        
        error('Not Yet Complete')
        
    otherwise
        
        error('err')
        
end

nRows_RowPartition_1 = nchoosek(m + n - k + 2, 2);
nRows_RowPartition_2 = nchoosek(m + o - k + 2, 2);

R1 = Sk(1 : nRows_RowPartition_1, :);
R2 = Sk(nRows_RowPartition_1 + 1 : end , : );

max_R1 = max(max(abs(R1(R1~=0))));
min_R1 = min(min(abs(R1(R1~=0))));

max_R2 = max(max(abs(R2(R2~=0))));
min_R2 = min(min(abs(R2(R2~=0))));

r1 = max_R1 ./ min_R1;
r2 = max_R2 ./ min_R2;

fprintf('Ratio r1 : %e \n',r1);
fprintf('Ratio r2 : %e \n', r2);

end



function Sk = BuildSylvesterMatrix_3Polys_3Eqns(fxy, gxy, hxy, m, n, o, k)
% BuildSylvesterMatrix_3Polys(fxy, fxy2, gxy, hxy, m, n, k)
%
% Build the Sylvester subresultant matrix S_{k}(f,g)
%
%
% Inputs
%
% fxy : (Matrix)
%
% fxy2 : (Matrix)
%
% gxy : (Matrix)
%
% hxy : (Matrix) Coefficients of the polynomial f(x,y), g(x,y) and
% h(x,y)
%
% m : (Int) Degree of f(x,y)
%
% n : (Int) Degree of g(x,y)
%
% o : (Int) degree of h(x,y)
%
% k : Index of the Sylvester subresultant matrix S_{k} to be constructed.


global SETTINGS
switch SETTINGS.SYLVESTER_MATRIX_VARIANT
    case 'T'
        
        
        T = BuildT_3Polys_3Eqns(fxy, gxy, hxy, m, n, o, k);
        
        Sk = T;
        
    case 'DT'
        
        D = BuildD_3Polys_3Eqns(m, n, o, k);
        T = BuildT_3Polys_3Eqns(fxy, gxy, hxy, m, n, o, k);
        
        % Build the matrix DT
        Sk = D*T;
        
    case 'DTQ'
        
        
        D = BuildD_3Polys_3Eqns(m, n, o, k);
        T = BuildT_3Polys_3Eqns(fxy, gxy, hxy, m, n, o, k);
        Q = BuildQ_3Polys(m, n, o, k);
        
        % Build the matrix DTQ
        Sk = D*T*Q;
        
        
    case 'TQ'
        
        T = BuildT_3Polys_3Eqns(fxy, gxy, hxy, m, n, o, k);
        Q = BuildQ_3Polys(m, n, o, k);
        
        
        Sk = T * Q;
        
    case 'DTQ Denominator Removed'
        
        error('Not Yet Complete')
        
    otherwise
        
        error('err')
        
end


end
