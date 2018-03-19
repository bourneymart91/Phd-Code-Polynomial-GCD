
function intercept = ComputeXIntercept(CP,degree)
% Calculate the intercept of the X axis of the line between the first and
% last control point.
%
% Inputs.
%
% CP :  Set of control points [x,y]
%
% degree :  Degree of control points.
%
%
%
%                       Outputs
%
% intercept :   Intercept of the x axis between C(0,0) and C(n,n)
%
%

% set first control point to be (x0,y0)
x_0 = CP(1, 1);
y_0 = CP(1, 2);

% set last control point to be (x1,y1)
x_1 = CP(degree + 1, 1);
y_1 = CP(degree + 1, 2);

% obtain change in y
d_y = y_1 - y_0;

% obtain change in x
d_x = x_1 - x_0;

% obtain gradient
m = d_y ./ d_x;

% calculate intercept of x axis by the line between (x0,y0) and (x1,y1)
intercept = x_1 - (y_1./m) ;

end