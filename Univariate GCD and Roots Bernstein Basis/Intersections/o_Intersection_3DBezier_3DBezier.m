function [] = o_Intersection_3DBezier_3DBezier(ex_num)
% Get the points of intersection of two Bezier curves defined in 3
% Dimensional space

% % Get the control points of the first Curve
% [x y z]
CP1 = Examples_3DBezier_ControlPoints(ex_num)

% Get the polynomials x(t), y(t), and z(t)
xt_1 = CP1(1,:);
yt_1 = CP1(2,:);
zt_1 = CP1(3,:);


% Add noise indepenedent of the value of the coefficients
max_error = 1e-1;
error_grid = rand(size(CP1)) .* max_error;

% Define the Control points of the secnod curve, to be the first curve with
% minor perturbations.
CP2 = CP1 + error_grid;

% Get the polynomials x(t), y(t), and z(t)
xt_2 = CP2(1,:);
yt_2 = CP2(2,:);
zt_2 = CP2(3,:);

% Get the number of control points.
[~,nControlPoints] = size(CP2)

% Get the degree of the curve.
m = nControlPoints - 1;

% Print out the curves.
fprintf('Control Points of Curve C_{1}: \n')
display(CP1)
fprintf('Control Points of Curve C_{2}: \n')
display(CP2)

% % Graph the curves

% For Graphing the curve
t = linspace (0, 1, 100);
Q3D_1 = Bezier (CP1, t);
Q3D_2 = Bezier (CP2, t);


figure('name','Bézier Curve Intersection in 3D')
% Plot Bézier curve
plot3(Q3D_1 (1, :), Q3D_1 (2, :), Q3D_1 (3, :), ' b'),
hold on
plot3(Q3D_2 (1, :), Q3D_2 (2, :), Q3D_2 (3, :), ' b'),
% plot control polygon
plot3( CP1(1, :), CP1 (2, :), CP1 (3, :), ' g : ')
plot3( CP2(1, :), CP2 (2, :), CP2 (3, :), ' g : ')
% plot control points
plot3 (CP1 (1, :), CP1 (2, :), CP1 (3, :), ' ro')
plot3 (CP2 (1, :), CP2 (2, :), CP2 (3, :), ' ro')
hold off

% % Implicitize the curve C2
x = sym('x')
y = sym('y')
z = sym('z')

C1_f = BuildC1_sym(sym(xt_1),m,x)
C1_g = BuildC1_sym(sym(yt_1),m,y)
C1_h = BuildC1_sym(sym(zt_1),m,z)

Sylvester = ...
    [
    C1_f            zeros(2*m,m+1)          C1_g
    zeros(2*m,m+1)    C1_f                  C1_h
    ];
det(Sylvester)
% Calculate their intersections

end

function C1 =  BuildC1_sym(f,m,x)
% Build the symbolic toeplitz matrix C1, containing coefficients of
% polynomial f.


% Get Binomials corresponding to the basis elements of f.
Bi_m = GetBinomials(m);


% Transpose f so that it is a row vector
f = f';

% Get f in the modified bernstein basis
f = (f - x) .* Bi_m;

% for each column in C1 j=0,...,m
for j = 1:1:m+1
    
    % Insert the coefficients of f(w) into the column.
    C1(j:m+j,j) = (f.*Bi_m);
end


end