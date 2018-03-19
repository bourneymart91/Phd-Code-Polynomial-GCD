function [] = o_Intersection_ImplicitAndBezier()
% Get the points of intersection between an implicitly defined curve C_{1} 
% and a Bézier curve C_{2}


% Get the coefficients of the implicitly defined curve
x = sym('x');
y = sym('y');

% Define the First curve (Implicit curve)
ex_num_implicit = '1';
switch ex_num_implicit
    case '1'
        C1 = 3*(x)-2*y.^2;
    case '2'
        
        C1 = y+x-1.3;
end
% Get the parametric equation of the bezier curve.
ex_num_Bezier = '1';
switch ex_num_Bezier
    case '1'
        C2_xt = [0 ;0.9; 0.7; 0.3];
        C2_yt = [1.1; 1.2; 0.6; 0.5];
    case '2'
        C2_xt = [0  ; 1/3; 2/3; 1];
        C2_yt = [1.1; 0.6; 0.8; 1.5]
    case '3'
        C2_xt = [0.1; 0.5];
        C2_yt = [0.8; 1.2];
        
end

% Get the interval over which the bezier curve is defined
min_t = min(C2_xt)
max_t = max(C2_yt)

% Get the matrix of control points
cp_f = [C2_xt' ; C2_yt'];

% Given a set of control points, plot the bezier curve.
degree_f = size(cp_f,2) - 1;
t = (min_t:0.01:max_t);

% Get a set of points for the Bezier Curve
pts_f = 0;
for i = 0:1:degree_f
    pts_f = pts_f + kron( nchoosek(degree_f,i).*((1-t).^(degree_f - i)) .* (t.^(i)) ,cp_f(:,i+1));
end

% Plot the implicit and the Bézier curve
figure()
hold on
%plot(t.^2,t.^2+5,'g')
ezplot(C1,[min_t,max_t])
plot(pts_f(1,:),pts_f(2,:),'b')
hold off

% Get the polynomial whose roots are to be found
poly = subs(C1,{x y},{C2_xt C2_yt});
poly = double(poly);

% Caculate roots by matlab method
vRoots = o_roots_matlab(poly);

% % 
% Define the global variables for my root finder.
global SETTINGS
SETTINGS.PLOT_GRAPHS = true;
SETTINGS.BOOL_LOG = false;
SETTINGS.BOOL_ALPHA_THETA = true;
SETTINGS.MEAN_METHOD = 'Geometric Mean Matlab Method';

% Set the sylvester matrix variant to be used. 
% 'T'
% 'DT'
% 'TQ'
% 'DTQ'
SETTINGS.SYLVESTER_MATRIX_VARIANT = 'DTQ';

SETTINGS.SNTLN_METHOD = 'None';
SETTINGS.APF_METHOD = 'None';

SETTINGS.APF_BUILD_METHOD = 'Standard APF' ;
SETTINGS.DECONVOLUTION_METHOD = 'single';

% Calculate roots by my method
%vRoots = o_roots_mymethod(poly);

% % Print out the roots and the corresponding intersection points.

[num_roots,~] = size(vRoots);

for i = 1:1:num_roots
    % Get the ith root
    r = vRoots(i);
    fprintf('Root : %2.4f \n', r)
    
    % Get the x coordinate
    x_ord = BernsteinEvaluate(C2_xt, r);
    % Get the y coordinate
    y_ord = BernsteinEvaluate(C2_yt, r);
    fprintf('Intersection Point : (%2.4f, %2.4f) \n',x_ord,y_ord)
    
end

end


