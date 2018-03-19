function [] = o_Intersection_3DBezierCurve_Surface
% Get Intersection Between a Bézier Curve in 3D space and a surface defined
% implicitly.

% A Bézier Curve in 3 Dimensions

% An Implicitly defined surface


%% The Curve defined by Control Points
CP = Examples_3DBezier_ControlPoints('2');

%
% For Graphing the curve
t = linspace (0, 1, 100);
Q3D = Bezier (CP, t);



%%
% Get the parametric functions x(t), y(t) and z(t) for the Bézier Curve
xt = CP(1,:);
yt = CP(2,:);
zt = CP(3,:);

% Initialise Symbolic Variables
x = sym('x');
y = sym('y');
z = sym('z');


%%
figure()

% Plot Bézier curve
plot3(Q3D (1, :), Q3D (2, :), Q3D (3, :), ' b'),
hold on

% plot control polygon
plot3(CP (1, :), CP (2, :), CP (3, :), ' g : ')

% plot control points
plot3 (CP (1, :), CP (2, :), CP (3, :), ' ro')

% Plot Implicit Surface
ezsurf(x,y,x+y);

xlim([-2,2])
ylim([-2,2])
zlim([-2,2])
view (3);
box;
%%

S2 = x + y - z;

S1_funx = x;
S1_funy = y;
S1_funz = x + y;


S3 = subs(S2,{x,y,z},{xt,yt,zt});

S3 = double(S3)';

% Get the set of roots by matlab method
vRoots = o_roots_matlab(S3);


%roots(flipud(S3))

% Get number of roots in vRoots
[nRoots,~] = size(vRoots);
for i = 1:1:nRoots
   % Get the root
   r = vRoots(i);
   % For each root, evaluate x(t), y(t) and z(t)
   xval = BernsteinEvaluate(xt',r);
   yval = BernsteinEvaluate(yt',r);
   zval = BernsteinEvaluate(zt',r);
   
   display([xval yval zval]);
end




end