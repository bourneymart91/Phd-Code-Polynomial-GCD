function [] =  o_Intersection_2DBezier_2DBezier(ex_num, emin, emax, ...
    mean_method, bool_alpha_theta, low_rank_approx_method, apf_method)
% o_Intersection_2DBezier_2DBezier(ex_num, emin, emax, mean_method, ...
%       bool_alpha_theta, low_rank_approx_method, apf_method)
%
% Given an example number, get the control points of 2 Bezier Curves, and
% calculate their intersections.
%
% Inputs.
% 
% ex_num :
%
% emin:
%
% emax:
%
% mean_method:
%
% bool_preproc :
%
% low_rank_approx_method : 
%
% apf_method :
%
% 
% % Example
%
% >> o_Intersection_2DBezier_2DBezier('1',1e-12,1e-11,'Geometric Mean
% Matlab Method', 'y','Standard STLN','None')


% 1. Get Control Points of polynomial f.
% 2. Get Control Points of polynomial g.
% 3. Implicitize Polynomial g.
% 4. Substitute the parametric expressions of f for x and y in t, into the
%    implicit polynomial g.

close all; clc; 


sylvester_matrix_variant = 'DTQ';
rank_revealing_metric = 'Minimum Singular Values';
deconvolution_method_hx = 'Batch';
deconvolution_method_wx = 'Batch';
deconvolution_preproc = true;


SetGlobalVariables_Roots(ex_num, emin, emax, ...
    mean_method, bool_alpha_theta, low_rank_approx_method, apf_method,...
    sylvester_matrix_variant, rank_revealing_metric, deconvolution_method_hx,...
    deconvolution_method_wx, deconvolution_preproc);

switch ex_num
    case '1'
        ex_num_f = '1';
        ex_num_g = '2';
end


% Get a set of control points for the planar Bezier curve f
CP_f = Examples_2DBezier_ControlPoints(ex_num_f);

% Get a set of control points for the planar Bezier curve g
CP_g = Examples_2DBezier_ControlPoints(ex_num_g);

array_cp{1} = CP_f;
array_cp{2} = CP_g;

Plot2DBezier(array_cp);

% %
% %
% %
% Implicitize the Bezier Control Points of g

fprintf('Implicit representation of g:\n')
[gxy, symbolic_expression] = ImplicitizeBezierBySylvester(CP_g);


figure()
hold on
xmin = 0;
xmax = 1.2;
ymin = -1;
ymax = 2;
ezplot(symbolic_expression , [xmin,xmax,ymin,ymax])
hold off

% Print the coefficients of the polynomial g(x,y) in Bernstein form.
PrintCoefficients_Bivariate_Bernstein(gxy,'C_{2}:')

% %
% %
% %
% Get the parametric expressions of polynomial f in the brn basis



% Put the sets of control points into an array
CP_arr{1} = CP_f;
CP_arr{2} = CP_g;

% Get coefficients of f(x)
f_x =  CP_f(1,:)';

% Get Coefficients of f(y)
f_y =  CP_f(2,:);

% Get Coefficients of the Bernstein polynomial obtained by susbstituting
% x(t) and y(t) from f(x,y) into g(x,y).


vCoefficients_Bernstein_Poly = Substitute(CP_f, gxy);
m = GetDegree(vCoefficients_Bernstein_Poly);
x_vec = 0:(1/(m)):1;

arrCP{1} = [x_vec; vCoefficients_Bernstein_Poly'];

Plot2DBezier(arrCP);


% % Get the roots of the polynomial
roots = o_roots_mymethod(vCoefficients_Bernstein_Poly);
%roots = o_roots_Matlab(coef_Bernstein_Poly);

% % 
% Given roots are calculated, plug in roots to one of the original 
% parametric equations
[nRoots,~] = size(roots);

% for each root
for i = 1:1:nRoots
    r_i = roots(i);
    if isreal(r_i)
        x = Bernstein_Evaluate(f_x,r_i);
        y = Bernstein_Evaluate(f_y',r_i);
        fprintf('The root %2.4e gives intersection point %2.3f , %2.3f \n',r_i,x,y)
   
    end
end
end



