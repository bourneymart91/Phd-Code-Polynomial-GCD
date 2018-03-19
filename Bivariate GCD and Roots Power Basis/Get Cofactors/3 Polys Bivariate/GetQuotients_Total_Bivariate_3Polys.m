function [uxy, vxy, wxy] = GetQuotients_Total_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, k)
% Given two polynomials and the knowledge of the degree of the GCD. Obtain
% the two quotient polynomials u and v.
%
% % Inputs.
%
% fxy : (Matrix) Coefficients of the polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of the polynomial g(x,y)
%
% hxy : (Matrix) Coefficients of the polynomail h(x,y)
%
% m n o : (Int) (Int) (Int) Total degree of f(x,y), g(x,y) and h(x,y)
%
% k : (Int) Total degree of GCD d(x,y)
%
% % Outputs
%
% uxy : (Matrix) Coefficients of polynomial u(x,y)
%
% vxy : (Matrix) Coefficients of polynomial v(x,y)
%
% wxy : (Matrix) Coefficients of polynomial w(x,y)


% %
% %
% Replace fxy with the padded version
% Get degree of f(x,y) with respect to x and y
[m1, m2] = GetDegree_Bivariate(fxy);
padd_matrix = zeros(m+1,m+1);
padd_matrix(1:m1+1,1:m2+1) = fxy;
fxy = padd_matrix;

% %
% %
% Replace gxy with the padded version
% Get degree of g(x,y) with respect to x and y
[n1, n2] = GetDegree_Bivariate(gxy);
padd_matrix = zeros(n+1,n+1);
padd_matrix(1:n1+1,1:n2+1) = gxy;
gxy = padd_matrix;


%
%
[o1, o2] = GetDegree_Bivariate(hxy);
padd_matrix = zeros(o+1,o+1);
padd_matrix(1:o1+1,1:o2+1) = hxy;
hxy = padd_matrix;


% Build the kth Sylvester Subresultant matrix S_{k}(f,g,h)
Sk = BuildT_Total_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, k);

% Get the optimal column for removal from S(f,g)
idx_optColumn = GetOptimalColumn_Total(Sk);


% % Having found the optimal column, obtain u and v the quotient polynomials.
Akj = Sk;
cki = Sk(:,idx_optColumn);
Akj(:,idx_optColumn) = [];

x_ls = SolveAx_b(Akj,cki);

% Obtain the solution vector x = [-v;u]
vecx =[
    x_ls(1:(idx_optColumn)-1);
    -1;
    x_ls(idx_optColumn:end);
    ]  ;

% Get number of coefficients in u(x,y) and v(x,y)
nCoefficients_vxy = nchoosek(n-k+2,2);
nCoefficients_wxy = nchoosek(o-k+2,2);
nCoefficients_uxy = nchoosek(m-k+2,2);

% get the vector of coefficients of v
vxy_calc = vecx(1:nCoefficients_vxy);
wxy_calc = vecx(nCoefficients_vxy + 1 : nCoefficients_vxy + nCoefficients_wxy);
uxy_calc = (-1) .* vecx(nCoefficients_vxy + nCoefficients_wxy + 1 : end);


% Get u(x,y) and v(x,y) in matrix form
try
    zeros_vww = zeros(nchoosek(n-k-1+2, 2), 1);
catch
    zeros_vww = 0;
end
try
    zeros_uww = zeros(nchoosek(m-k-1+2, 2), 1);
catch
    zeros_uww = 0;
    
end
try
    zeros_www = zeros(nchoosek(o-k-1+2, 2), 1);
catch
    zeros_www = 0;
end


% Get coefficients of u(x,y), v(x,y) and w(x,y) as matrices
uxy = GetAsMatrix_Version1([uxy_calc; zeros_uww], m-k, m-k);
wxy = GetAsMatrix_Version1([wxy_calc; zeros_www], o-k, o-k);
vxy = GetAsMatrix_Version1([vxy_calc; zeros_vww], n-k, n-k);





end
