function [ux, vx] = GetCofactorsCoefficients_2Polys(fx, gx, k)
% Get Quotient polynomials u(x) and v(x) such that
% f(x)/u(x) = g(x)/v(x) = d(x)
%
% % Inputs.
%
% fx  : Coefficients of polynomial f(x)
%
% gx  : Coefficients of polynomial g(x)
%
% k     : Degree of common divisor d(x) of degree k.
%
% % Outputs
%
% ux : Coefficients of the polynomial u(x)
%
% vx : Coefficients of the polynomial v(x)


% Get degree of polynomial g(x).
n = GetDegree(gx);

% Build the t-th subresultant S_{t}(f,g)
Sk = BuildT(fx, gx, k);

% Get index of optimal column for removal from S_{k}(f,g)
[~, idx_col] = GetMinDistance(Sk);

% Get the matrix A_{t}(f,g), S_{t} with the optimal column removed
Ak = Sk;
Ak(:,idx_col) = [];

% Get the column c_{t} removed from S_{t}
ck = Sk(:,idx_col);

% get the vector x_ls
x_ls = SolveAx_b(Ak,ck);

% insert 1 into x_ls in the position corresponding to the col
x = [x_ls(1:idx_col-1); -1 ; x_ls(idx_col:end)];

% Split x into v(x) and u(x)

% Get the number of coefficients in v(x)
nCoeffs_vx = n-k+1;

% Get coefficients of v(x) and u(x)
vx =  x(1:nCoeffs_vx);
ux = -x(nCoeffs_vx + 1:end);


end

