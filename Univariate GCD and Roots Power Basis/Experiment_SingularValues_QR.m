
function [] = Experiment_SingularValues_QR
% Experiment on the Singular values and QR Decomposition of a matrix

% Get the roots of two polynomials f(x) and g(x)
roots_fx = Examples_Roots('1');
roots_gx = Examples_Roots('1b');

% Print the factorization of the two polynomials f(x) and g(x)
PrintFactorization(roots_fx,'f')
PrintFactorization(roots_gx,'g')

% Get the coefficients of the two polynomials f(x) and g(x)
fx = GetCoefficients(roots_fx)
gx = GetCoefficients(roots_gx)

% Add noise to the coefficients of one of the polynonmials
el = 1e+1;
gx = gx + el.*rand(size(gx))

% Get the degree of polynomial f(x)
[nRows_f,~] = size(fx);
m = nRows_f -1;

% Get the degree of polynomial g(x)
[nRows_g,~] = size(gx);
n = nRows_g -1;



C1 = BuildC1(fx,n,1);
C2 = BuildC1(gx,m,1);

S = [C1 C2]

% Perform SVD
[U,e,V] = svd(S)

% Perform QR 
[Q,R] = qr(S)

vDiagR = sort(abs(diag(R)),'descend');
vSingularValues = sort(diag(e),'descend');
vLowerDiag = sort(abs(diag(L)),'descend')

vDiagR = log10(vDiagR)
vSingularValues = log10(vSingularValues)
vLowerDiag = log10(vLowerDiag)

figure('name','Singular Values')
plot((vSingularValues),'b-s');
hold off

figure('name','Diagonal R Values')
plot((vDiagR),'b-s')
hold off
 
figure('name','LU')
plot(vLowerDiag,'b-s')
hold off

end