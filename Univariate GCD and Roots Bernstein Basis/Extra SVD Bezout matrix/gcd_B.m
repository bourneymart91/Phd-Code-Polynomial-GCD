
function [dx_QR, dx_lu] = gcd_B(Bez_fg, t, theta)

% This function computes the coefficients of the GCD of the polynomials
% f and g whose Bezout matrix is S.

% deg       :  The degree of the GCD of f and g.

% theta     :  The value of theta. Note that theta=1 for the Bernstein
%              basis, and theta~=1 for the modified Bernstein basis.

% CGCD_qr   :  The coefficients of the GCD computed from the QR 
%              decomposition of S.

% CGCD_lu   :  The coefficients of the GCD computed from the LU 
%              decomposition of S.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use the QR and LU decompositions to reduce S to upper triangle form.
[~,R] = qr(Bez_fg);

[~,U] = lu(Bez_fg);

% Calculate the coefficients of the GCD from the last nonzero row of R.
dx_QR = calculateCoefficients(R, t, theta);

% Calculate the coefficients of the GCD from the last nonzero row of U.
dx_lu = calculateCoefficients(U, t, theta);