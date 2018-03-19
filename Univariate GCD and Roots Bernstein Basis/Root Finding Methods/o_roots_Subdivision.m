function [ root_mult_matrix ] = o_roots_Subdivision(fx)
% Given the coefficients of the polynomial in Bernstein form, compute the
% roots using the subdivision method by Schneider in Graphics Gems.
%
% Inputs.
%
% fx : (Vector) Coefficients of polynomial f(x)
%
% Outputs.
%
% root_mult_matrix : (Matrix) Contains roots and corresponding multiplicity
%

% Global variables
global MAXDEPTH
MAXDEPTH = 30;

global EPSILON
EPSILON = 1e-25;

% Initialise iteration number 
ite_num = 1;

% Get degree of polynomial f(x)
m = GetDegree(fx);

% Get the set of control points of the curve
a = 0;
b = 1;

% Initialise a matrix of control points, where each row contains the (x,y)
% pair.
Pk = zeros(m + 1, 2);

% Get the control points
for i = 0 : 1 : m
    % Get the x ordinate
    x_ord = a + (i/m).*(b-a);
    Pk(i+1,:) = [x_ord fx(i+1)];
end

% Plot the coefficients and the control points of the curve
figure_name = sprintf('Bernstein Polynomial Curve %i',ite_num - 1);
Plot_fx(fx, a, b, figure_name);

% Set the degree of polynomial f
degree = m;

W_DEGREE = m;

% Set the intial depth to 1
ini_depth = 1;

% Set the number of found roots
root_count = 0;

% set the current reached depth to 1
reached_depth = 1;

% Obtain the roots using the bisection method.
[t,~,reached_depth] = FindRoots(Pk, degree, ini_depth, W_DEGREE, root_count, reached_depth);

root_mult_matrix = [];

if isempty(t)
    fprintf('\n')
    fprintf('ROOTS CALCULATED BY SUBDIVISION FUNCTION \n');
    fprintf('No Roots were Found\n')
    return
end

% %
% Given the set of calculated simple roots produce a polynomial in the 
% bernstein form.
new_roots = [real(t(:,1)), ones(length(t),1)];
root_mult_matrix = new_roots;

% Get coefficients of polynomial r1 in the scaled bernstein form
r1_bi = BuildPolyFromRoots(new_roots);

% Divide by binomial coefficients to obtain r1 in standard bernstein basis.
r1_x = GetWithoutBinomials(r1_bi);

% Obtain the polynomial f_{2}, given by the removal of the calculated roots
% from f_{1}.
fx = Deconvolve(fx, r1_x);

% Plot the coefficients and the control points of the curve
figure_name = sprintf('Bernstein Polynomial Curve %i',ite_num );
Plot_fx(fx, a, b, figure_name);

fprintf('\n\n')


% While f_{2} is not a scalar.
while GetDegree(fx) >= 1
    
    % Increment iteration number
    ite_num = ite_num + 1;
    
    % Get degree of polynomial f2
    m_ite = GetDegree(fx);
    
    % Get the control points of f2
    Pk = GetControlPoints(a, b, fx);
    
    % Plot f2
    figure_name = sprintf('Bernstein Polynomial Curve %i',  ite_num);
    Plot_fx(fx, a, b, figure_name);
    
    % Get roots of polynomial f(x)
    [t_new,~,reached_depth] = FindRoots(Pk, m_ite, 1, m_ite, 1, 1);
    
    % add the found roots to a list of roots
    root_mult_matrix = [root_mult_matrix; [t_new ones(length(t_new),1)] ];
    
    % if no new roots are found, end this while loop
    [nEntries_t,~] = size(t_new);
    if (nEntries_t == 0)
        break;
    end
    
    % Build Polynomial from newly obtained roots
    new_roots = [real(t_new(:,1)), ones(length(t_new),1)];

    % Get polynomial coefficients of r1 in scaled bernstien basis
    r1_bi = BuildPolyFromRoots(new_roots);
    
    % Get r1 in standard bernstein basis.
    r1_x = GetWithoutBinomials(r1_bi);
    
    % Perform deconvolution and obtain remaining part of original input
    % polynomial. Update f2 to be the result of the deconvolution.
    fx = Deconvolve(fx, r1_x);
    
    
    
    
end

end









