function [fxy, gxy, uxy, vxy, dxy, m, n, n_t, m_t, t, t1, t2] = GetCoeff(ex)
% Given the roots of the polynomials f(x,y) g(x,y), u(x,y), v(x,y) and
% d(x,y) return the matrices of their coefficients. The entries of these
% matrices are of the form a_{i,j} x^{m1-i} y^{m2-j}. The first entry has 
% the maximum powers x^{m1}y^{m2}, and the last bottom right entry has 
% powers x^{0}y^{0}.

% Get the polynomial roots and multiplicities
[roots_f_x, roots_f_y, ...
    roots_g_x, roots_g_y, ...
    roots_u_x, roots_u_y, ...
    roots_v_x, roots_v_y, ...
    roots_d_x, roots_d_y,...
    m, n, m_t, n_t,...
    t,t1,t2] = Example_SeparableRoots(ex);

% % BUILD THE POLYNOMIALS
% Get coefficients of polynomial f 
f_x_poly_pwr = BuildPoly_Pwr(roots_f_x);
f_y_poly_pwr = BuildPoly_Pwr(roots_f_y);
fxy = f_x_poly_pwr * f_y_poly_pwr';

% Get coefficients of polynomial g
g_x_poly_pwr = BuildPoly_Pwr(roots_g_x);
g_y_poly_pwr = BuildPoly_Pwr(roots_g_y);
gxy = g_x_poly_pwr * g_y_poly_pwr';

% Get coefficients of polynomial u
u_x_poly_pwr = BuildPoly_Pwr(roots_u_x);
u_y_poly_pwr = BuildPoly_Pwr(roots_u_y);
uxy = u_x_poly_pwr * u_y_poly_pwr';

% Get coefficients of polynomial u
v_x_poly_pwr = BuildPoly_Pwr(roots_v_x);
v_y_poly_pwr = BuildPoly_Pwr(roots_v_y);
vxy = v_x_poly_pwr * v_y_poly_pwr';

% Get Coefficients of polynomial d
d_x_poly_pwr = BuildPoly_Pwr(roots_d_x);
d_y_poly_pwr = BuildPoly_Pwr(roots_d_y);
dxy = d_x_poly_pwr * d_y_poly_pwr';

end

