function [] = o_Bezout(pz,qz)
% Given the two polynomials fx and gx in bernstein form, compute the degree
% of the GCD by the bezout matrix method.


% % This example matches that in the paper "Bini Bernstein-Bezoutian matrices". 
% % Given coefficients in power basis, conver to bernstein basis
%   pz = [4 0 -5 0 1];
%   qz = [1/2 -1/4 -2 1];
% % 
%   pz = PowerToBernstein(pz);
%   qz = PowerToBernstein(qz);


% Get the degree of each of the polynomials
m = length(pz)-1;
n = length(qz)-1;

% Get the maximum degree of the polynomials
max_deg = max(m,n);

% Degree Elevate the Polynomial coefficients to the maximum degree
pz = DegreeElevate(pz,max_deg-m);
qz = DegreeElevate(qz,max_deg-n);

% Build the bezoutian matrix
B = Bezout(pz,qz);

% perform svd on the bezoutian matrix

vec_SingVal = svd(B);

% Plot the singular values
figure()
hold on
title('Singular Value Decomposition of the Bezoutian Matrix')
plot(log10(vec_SingVal),'-s')
hold off


% Get the QR decomposition of the bezoutian
[q,r] = qr(B);

% Since B is always square, then R is upper triangular
% Get the diagonal entries of R
diag_r = diag(r);

% plot the diagonal entries
figure()
hold on
title('Diagonal entries of R from the QR decomposition of the Bezoutian Matrix')
plot(log10(abs(diag_r)),'-s')
hold off


% Further work, as seen in paper "Bini Bernstein-Bezoutian matrices"

% B = [B_k P ; Q R]
% get the principle subresultant B_{2} of B, which is the 2x2 submatrix
k = 3;
B_k = B(1:k,1:k);
% Let P be the remaining n-j columns of the first j rows
P = B(1:k,k+1:end);
% Let Q be the remaining n-j rows under B_k
Q = B(k+1:end, 1:k);
% Let R be the remaining n-j columns n-j rows
R = B(k+1:end, k+1:end);

Schur = R - Q*pinv(B_k)*P;


