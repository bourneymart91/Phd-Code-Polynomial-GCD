function [] = o_IntersectionImplicitParametric()
% Get the intersection between an Implicitly defined surface S_{1}(x,y,z) and a
% parametrically defined surface.

% S_Explicit
% An Explicitly defined surface z = f(x,y)


% Initialise the symbolic parametric variables s and t
s = sym('s');
t = sym('t');

% Initialise the symbolic variables x y and z
x = sym('x');
y = sym('y');
z = sym('z');


% The surface S_{1} is implicitly given by S_{1} = x + y - z.
% f(x,y,z) = 0
S1 = x + y - z;

% The parametric definitions of S_{1} are given as follows, and are used in
% plotting.
S1_fun_x = s;
S1_fun_y = t;
S1_fun_z = s + t;

% Get the parametrically defined surface S2
S2_fun_x_coef = ...
    [
        5   0
        1   0
    ];
S2_fun_x = GetSymBivariatePoly(S2_fun_x_coef);

S2_fun_y_coef = ...
    [
        0   1
        0   0
    ];
S2_fun_y = GetSymBivariatePoly(S2_fun_y_coef);

S2_fun_z_coeff = ...
    [
        0   1
        1   1
    ];
S2_fun_z = GetSymBivariatePoly(S2_fun_z_coeff);

% Plot the surfaces
figure('name','Surface Plot')
hold on
Surf1 = ezsurf(S1_fun_x,S1_fun_y,S1_fun_z);
Surf2 = ezsurf(S2_fun_x,S2_fun_y,S2_fun_z);
hold off

% Substitute x(s,t),y(s,t) and z(s,t) from S2 into S1(x,y,z) to obtain a
% bivariate polynomial h(s,t).
h_st = subs(S1,{x y z},{S2_fun_x_coef, S2_fun_y_coef,S2_fun_z_coeff});
h_st = double(h_st);

% Get degree of h(s,t)
[m1, m2] = GetDegree_Bivariate(h_st);
m = m1 + m2;

% Get the roots of h(s,t).
o_roots_mymethod(h_st,m)



end
