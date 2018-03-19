function [C1_xy_matrix] = Implicitize_Rational_Parametric_BySylvesterMatrix ...
    (C1_xt,C1_yt,C1_wt)

% Implicitize the rational parametric curve C1 = x(t)/w(t),y(t)/w(t).
% Return the matrix of coefficients of C1(x,y)

% Get degree of x(t)
[r,~] = size(C1_xt);
m = r-1;

% Get degree of y(t)
[r,~] = size(C1_yt);
n = r-1;

% Initialise the symbol x
x = sym('x');
y = sym('y');

rhs_vec_x = x.*C1_wt;
rhs_vec_y = y.* C1_wt;

% Set the first coefficient of x(t) (corresponding to the constant term of the 
% polynomial) to a_{0} - x
f    = sym(C1_xt);
f = f - rhs_vec_x


% Set the first coefficient of y(t) to b_{0}-y
g    = sym(C1_yt);
g    = g - rhs_vec_y;


% Build a zero matrix C1
T1 = zeros(m+n,n);
T1 = sym(T1);

% Build the matrix C1
for i = 1:1:n
    T1(i:i+n,i) = f;
end

% Build a zero matrix C2
T2 = zeros(n+m,m);
T2 = sym(T2);

% Build the matrix C2
for i = 1:1:m
    T2(i:i+n,i) = g;
end


% Build the matrix C2
S = [T1 T2];

%fprintf('The implicit representation is given as');
symbolic_det_expression = expand(det(S))

% Split the poly into a series of polynomials in terms of y
cx_vec = coeffs(symbolic_det_expression, x);

% a degree 3 polynomial has a bivariate implicit representation of degree
% (3,3) = (t1,t2)

max_x = n;

% get highest power of y
max_y = n;

mat = zeros(max_x+1,max_y+1);

% for each column in cx_vec
[~,c] = size(cx_vec);

for i = 1:1:c
    
    % Get the coefficients of the polynomial in y, corresponding to x^{i}
    % x^{i} * (y + y^{2} + ... )
    
    coef_vec_y = fliplr(sym2poly(cx_vec(1,i)));
    
    [~,cy] = size(coef_vec_y);
    
    mat(i,1:cy) = coef_vec_y;
    
    
end

C1_xy_matrix = mat;
