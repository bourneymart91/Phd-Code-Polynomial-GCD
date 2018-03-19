function DYQ = BuildDYQ_STLN(x,m,n,k)
% BuildY(x,m,n,k)
%
% Build the matrix Y_{k}(x1,x2)
%
% % Inputs
%
% x : (Vector) Vector x = [x1 ; x2] where x1 and x2 are vectors.
%
% m : (Int) Degree of polynomial f(x,y)
%
% n : (Int) Degree of polynomial g(x,y)
%
% k : (Int) index of kth Sylvester subresultant.
%
% % Outputs
%
% DYQ : (Matrix) Matrix DYQ

% %
% Split the vector x into \hat{x}_{1} and \hat{x}_{2}

% Get number of coefficients in x1
nCoefficients_x1 = nchoosek(n-k+2,2);

% Get the number of zeros in a matrix containing entries of x1
nZeros_x1 = nchoosek(n-k+1, 2);

% Get the number of zeros in a matrix containing entries of x2
nZeros_x2 = nchoosek(m-k+1, 2);

% Split x into x1 and x2
x1 = x(1:nCoefficients_x1);
x2 = x(nCoefficients_x1 + 1 : end);

% Get vectors of coefficients of x_{v} x_{u} x_{1} and x_{2}
vec_x1 = [ x1 ; zeros(nZeros_x1, 1)];
vec_x2 = [ x2 ; zeros(nZeros_x2, 1)];

% Get the vectors as matrices of coefficients.
mat_x1 = GetAsMatrix(vec_x1, n-k, n-k);
mat_x2 = GetAsMatrix(vec_x2, m-k, m-k);

% Build the matrix T_{m}(\hat{x}_{1})
T1_x1 = BuildT1(mat_x1, n-k, m);
T1_x2 = BuildT1(mat_x2, m-k, n);

D = BuildD_2Polys(m, n-k);
Qm = BuildQ1(m);
Qn = BuildQ1(n);

DYQ = D * [T1_x1*Qm T1_x2*Qn] ;

end
