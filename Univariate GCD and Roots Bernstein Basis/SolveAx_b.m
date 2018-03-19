function x = SolveAx_b(A,b)
% Given the matrix A and the right hand side vector b, solve the problem
% Ax=b for x.
%
% Inputs.
%
% A : Coefficient matrix A
%
% b : Right hand side vector b

% %
% % Get x_ls by QR Decomposition method
% %
% %
% %

%warning('off')

[~,n2] = size(A);
[Q,R] = qr(A);
R1 = R(1:n2,:);
cd = Q'*b;
c = cd(1:n2,:);

x_ls_QR = R1\c;

% Calculate the residual
res_QR = b - (A*x_ls_QR);

% %
% %
% % Get x_ls by matlab pinv method
% %
%x_ls_SVD = pinv(A)*b;

%res_SVD = b - (A*x_ls_SVD);

% %

%if (res_QR < res_SVD)
    x = x_ls_QR;
%else
%    x = x_ls_SVD;
%end

%warning('on')

end