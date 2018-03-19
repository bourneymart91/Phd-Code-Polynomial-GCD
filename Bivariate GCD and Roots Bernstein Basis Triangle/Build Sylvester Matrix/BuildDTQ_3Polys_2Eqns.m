function DTQ = BuildDTQ_3Polys_2Eqns(fxy, gxy, hxy, m, n, o, k)
% BuildSubresultant(fxy_matrix_n,gxy_matrix_n,k1,k2,alpha,th1,th2)
%
% Build the sylvester subresultant matrix S_{k1,k2}.
%
% % Inputs.
%
%
% fxy1 : (Matrix) Coefficients of f(x,y)
%
% gxy : (Matrix) Coefficients of g(x,y)
%
% hxy : (Matrix) Coefficients of h(x,y)
%
% m : (Int)  : Degrees of polynomials f(x,y), g(x,y) and h(x,y)
%
% n : (Int) Degree of g(x,y)
%
% o : (Int) Degree of h(x,y)
%
% % Outputs
% 
% DTQ : Matrix DTQ


% Build the matrix D_{m1+n1-k1,m2+n2-k2} * T_{n1-k1,n2-k2}(f) * Q1_{n1-k1,n2-k2}
DT1Q1 = BuildDT1Q1(fxy, m, n - k);
DT2Q2 = BuildDT1Q1(fxy, m, o - k);


% Build the matrix D_{m1+n1-k1,m2+n2-k2} * T_{m1-k1,m2-k2}(g) * Q1_{m1-k1,m2-k2}
DT3Q3 = BuildDT1Q1(gxy, n, m - k);
DT4Q4 = BuildDT1Q1(hxy, o, m - k);


% Build the matrix DTQ
diagonal = blkdiag(DT1Q1,DT2Q2);
column = [DT3Q3; DT4Q4];

DTQ = [diagonal column];


end