function [uxy, vxy] = GetCofactors_2Polys(fxy, gxy, t)
% Given polynomials f(x,y) and g(x,y), and the degree of their GCD d(x,y).
% Compute the cofactors u(x,y) and v(x,y) by solving the Ax=b problem,
% where A is a modified Sylvester matrix with column removed, b is the
% removed column, and x is the vector of coefficients of u(x,y) and v(x,y)
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)
%
% t : (Int) Degree of polynomial d(x,y)
%
%
% % Outputs.
%
% uxy : (Matrix) Coefficients of polynomial u(x,y)
% 
% vxy : (Matrix) Coefficients of polynomial v(x,y)

% Get the degree of f(x,y) (Note m1 = m2)
[m1,~] = GetDegree_Bivariate(fxy);
m = m1;

% Get the degree of g(x,y) (Note n1 = n2)
[n1,~] = GetDegree_Bivariate(gxy);
n = n1;

% % Solve the Ax = b problem 

% Build the modified Sylvester matrix S_{k} = D^{-1}_{m+n-k}T_{k}(f,g)Q_{n-k,m-k}

Sk = BuildSylvesterMatrix_2Polys(fxy, gxy, m, n, t);

% %
% % Remove optimal column of S_{k}(f,g)

% Get index of optimal column for removal
opt_col = GetOptimalColumn(Sk);

% Get A_{k}, the matrix S_{k} with optimal column removed
Ak = Sk;
Ak(:,opt_col) = [];

% Get c_{k}, the column removed from S_{k}
ck = Sk(:,opt_col);

% Solve Ax = b
x_ls = SolveAx_b(Ak,ck);

% Insert '-1' into the opt_col position of x_ls to obtain x containing all
% coefficients of u(x,y) and v(x,y)
x = ...
    [...
        x_ls(1:opt_col-1) 
        -1
        x_ls(opt_col:end)
    ];

% %
% %
% Split x to obtain v(x,y) and u(x,y) and get u(x,y) and v(x,y) as matrices
% of coefficients


% If Q is included in Sylvester Matrix, then u(x,y) and v(x,y) coefficients
% are given by the vector x, otherwise, the vector gives coefficients in 
% scaled Bernstein form and must remove trinomial coefficients.
global SETTINGS
switch SETTINGS.SYLVESTER_MATRIX_VARIANT
    
    case 'T'
        Q = BuildQ_2Polys(m, n, t);
        x = Q\x;
    
    case 'DT'
        % Matrix Q is included in the solution vector so divide by this to
        % get coefficients u(x,y) and v(x,y)
        Q = BuildQ_2Polys(m, n, t);
        x = Q\x;
    
    case 'DTQ'
        % Q is included in S_{k} so u(x,y) and v(x,y) are 
    
    case 'TQ'
  
        
    case 'DTQ Denominator Removed'
        
        
    otherwise 
        error('SYLVESTER_MATRIX_VARIANT must be one of the following : \n T \n DT \n TQ \n DTQ')
end


% Get number of coefficients in v(x,y)
nCoefficients_vxy = nchoosek(n-t+2, 2);



% Get vectors of coefficients of u(x,y) and v(x,y)
v_vxy = x(1 : nCoefficients_vxy);
v_uxy = -1 .* x(nCoefficients_vxy+1:end);


% Get number of zeros in u(x,y)
nZeros_uxy = nchoosek(m-t+1,2);

% Get number of zeros in v(x,y)
nZeros_vxy = nchoosek(n-t+1,2);

% Get matrices of coefficients of u(x,y) and v(x,y)
uxy = GetAsMatrix([v_uxy ; zeros(nZeros_uxy, 1)], m - t, m - t);
vxy = GetAsMatrix([v_vxy ; zeros(nZeros_vxy, 1)], n - t, n - t);
    



end
