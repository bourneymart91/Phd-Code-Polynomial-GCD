function [] = o_Intersection_Implicit_Curves(bool_preproc, low_rank_approx_method)
% Given two implicit curves, calculate the points of intersection


% Get implicit curves C1 = f(x,y) = 0
% Get implicit curves C2 = g(x,y) = 0

SetGlobalVariables(bool_preproc,low_rank_approx_method);

% C1 := f(x,y) = 0.
% C2 := g(x,y) = 0.

% Define the symbolic variables x and y
x = sym('x');
y = sym('y');

% Get the curve f(x,y) = 0
curve1_Coefficients = ...
    [
    0 0 -1;
    0 0 0;
    1 0 0;
    ];
m = 2;

% Get the curve g(x,y) = 0
curve2_Coefficients = ...
    [
        0   -1
        1   0
    ];
n = 1;

% Get the symbolic representation of the curve f(x,y) = 0.
curve1_symbolic = GetSymbolicFromCoefficients(curve1_Coefficients,'x','y');
curve2_symbolic = GetSymbolicFromCoefficients(curve2_Coefficients,'x','y');



% Plot the curve h(x,y) = f(x,y)-g(x,y) = 0
curve3_symbolic = curve1_symbolic - curve2_symbolic;

% % 



% %
% % Plot the curves
figure('name','Plot')
hold on
ezplot(curve1_symbolic)
ezplot(curve2_symbolic)
ezplot(curve3_symbolic)
hold off


% %
[dxy] = o1(curve1_Coefficients,curve2_Coefficients,m,n)

% %
% Get the coefficients of the curve C3
C3 = ...
    [
    -5  1   1;
    -1  0   0;
    1   0   0;
    ];
% Get the degree of the curve C3
m = 2;


o_roots_mymethod(C3,m)



end