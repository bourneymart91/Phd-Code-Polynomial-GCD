
function boolFlatEnough = ControlPolygonFlatEnough(CP, degree)
% Given a control polygon, return a boolean as to whether it is flat enough
% to be considered a straight line.
%
% Inputs
%
%
% CP :  Set of control points which form the control polygon [x,y]
%
% degree :  Degree of control polygon
%
% Outputs
%
% boolFlatEnough : (Boolean)
%   1 : Control polygon is flat enough to be considered a line
%   0 : Control polygon is not flat enough.


global EPSILON


% Derive implicit equation for straight line connecting first and last 
% control point
a = CP(1,1) - CP(degree,1);

b = CP(1,2) - CP(degree,1);

c = (CP(1,1) * CP(degree,2)) - (CP(degree,1) * CP(1,2));

abSquared = (a.*a) + (b.*b);

for i = 1:1:degree
    
    distance(i) = a.*CP(i,1) + b.*CP(i,2) + c;
    
    if (distance > 0)
        distance(i) = (distance(i).*distance(i)) ./ abSquared;
    
    elseif (distance < 0)
        distance(i) = -((distance(i).* distance(i)))./ abSquared;
    end
end

max_distance_above = 0;
max_distance_below = 0;

for i =1:1:degree-1
    if (distance(i) < 0)
        fprintf('')
        max_distance_below = min(max_distance_below,distance(i));
    elseif (distance(i) > 0)
        max_distance_above = max(max_distance_above,distance(i));
    end
end

% Implicit equation for zero line
a1 = 0.0;
b1 = 1.0;
c1 = 0.0;

% implicit equation for "above" line
a2 = a;
b2 = b;
c2 = c+ max_distance_above;

det = (a1 .* b2) - (a2 .* b1);
dInv = 1.0 ./ det;

intercept_1 = (b1 * c2 - b2 * c1) * dInv;

% Implicit equation for "below" line

a2 = a;
b2 = b;
c2 = c+ max_distance_below;

det = (a1 .* b2) - (a2 .* b1);
dInv = 1.0 ./ det;

intercept_2 = (b1 * c2 - b2 * c1) * dInv;


% compute Intercepts of bounding box
left_intercept = min(intercept_1,intercept_2);
right_intercept = max(intercept_1,intercept_2);

error = 0.5 * (right_intercept-left_intercept);

if (error < EPSILON)
    
    %fprintf('Control Polygon is flat enough. \n')
    boolFlatEnough = true;
    
else
    %fprintf('Control Polygon is not flat enough. \n')
    boolFlatEnough = false;
end
end