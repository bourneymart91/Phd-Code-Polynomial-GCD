function [uxy, vxy, wxy] = GetQuotients_Both_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, k, k1, k2)
% GetQuotients(fxy_matrix,gxy_matrix,t1,t2,opt_alpha,th1,th2)
%
% Given two polynomials and the knowledge of the degree of the GCD. Obtain
% the two quotient polynomials u and v.
%
% % Inputs.
%
% fxy : Coefficients of the polynomial f(x,y)
%
% gxy : Coefficients of the polynomial g(x,y)
%
% hxy : Coefficients of the polynomial h(x,y)
%
% m n o : (Int) (Int) (Int) Total degrees of polynomials f, g and h
%
% k : (Int) Total degree of d(x,y)
%
% k1, k2 : (Int) (Int) : Degree of GCD d(x,y) with respect to x and y
%
% % Outputs
%
% uxy : (Matrix) Coefficients of polynomail u(x,y)
%
% vxy : (Matrix) Coefficients of polynomail v(x,y)
%
% wxy : (Matrix) Coefficients of polynomail w(x,y)


% Get the degree of polynomial f(x,y), g(x,y) and h(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);
[n1, n2] = GetDegree_Bivariate(gxy);
[o1, o2] = GetDegree_Bivariate(hxy);

% Build the Sylvester matrix S_{k,k1,k2}
Skk1k2 = BuildT_Both_Bivariate_3Polys(fxy, gxy, hxy, m, n, o, k, k1, k2);

% % Get the optimal column for removal.
idx_optColumn = GetOptimalColumn_Both(Skk1k2);

% Having found the optimal column, obtain u and v the quotient polynomials.
Atj = Skk1k2;
cki = Skk1k2(:,idx_optColumn);
Atj(:,idx_optColumn) = [];

% Get the solution vector.
x_ls = SolveAx_b(Atj,cki);

% Get number of non-zero entries in u(x,y) and v(x,y)
nNoneZeros_uxy = GetNumNonZeros(m1-k1, m2-k2, m-k);
nNoneZeros_vxy = GetNumNonZeros(n1-k1, n2-k2, n-k);
nNoneZeros_wxy = GetNumNonZeros(o1-k1, o2-k2, o-k);

% Get number of zero entries in u(x,y) and v(x,y)
nZeros_uxy = (m1-k1+1) * (m2-k2+1) - nNoneZeros_uxy;
nZeros_vxy = (n1-k1+1) * (n2-k2+1) - nNoneZeros_vxy;
nZeros_wxy = (o1-k1+1) * (o2-k2+1) - nNoneZeros_wxy;

% Obtain the solution vector x = [-v;u]
vecx =[
    x_ls(1:(idx_optColumn)-1);
    -1;
    x_ls(idx_optColumn:end);
    ]  ;

% Get the vector of coefficients of v(x,y)
vxy_calc = [...
            vecx(1:(nNoneZeros_vxy));
            zeros(nZeros_vxy,1)
          ];

% Get the vector of coefficients of w(x,y)
wxy_calc = [...
    vecx( nNoneZeros_vxy + 1 : nNoneZeros_vxy + nNoneZeros_wxy);
    zeros(nZeros_wxy,1)
    ];
      
% Get the vector of coefficients of u(x,y)
uxy_calc = [...
            (-1).*vecx( nNoneZeros_vxy + nNoneZeros_wxy + 1 :end);
            zeros(nZeros_uxy,1);
            ];
        

% Arrange coefficients of u(x,y), v(x,y) and w(x,y) as a matrix.
uxy = GetAsMatrix_Version1(uxy_calc, m1-k1, m2-k2);
wxy = GetAsMatrix_Version1(wxy_calc, o1-k1, o2-k2);
vxy = GetAsMatrix_Version1(vxy_calc, n1-k1, n2-k2);





end
