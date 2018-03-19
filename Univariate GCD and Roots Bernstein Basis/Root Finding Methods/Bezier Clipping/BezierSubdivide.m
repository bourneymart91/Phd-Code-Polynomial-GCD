
function [leftControlPoints, rightControlPoints] = BezierSubdivide(Pk, degree, c)
% Given a set of control points, subdivide such that two new sets of
% control points are obtained.
%
% % Inputs.
%
% Pk : (Matrix) Set of original control points.
%
% degree : (Int)  degree of control points
%
% c : (Float) Point of subdivision.
%
% % Outputs.
%
%
% leftControlPoints : (Matrix) The set of control points on the left of c.
%
% rightControlPoints : (Matrix) The set of control points on the right of c.
%

% obtain each set of control points
Pk_array{1} = Pk;

for i = 2 : 1 : degree + 1
    
    Pk_array{i} = deCasteljau(c, Pk_array{i-1})  ;
    
end

% From the generated sets of control points obtain set of control points
% for sub-curve P[a,c]
leftControlPoints = zeros(degree + 1, 2);

for i = 1 : 1 : degree + 1
    
    leftControlPoints(i,:) = Pk_array{i}(1,:);
    
end

% obtain set of control points for sub-curve P[c,b]
rightControlPoints = zeros(degree+1,2);
for i = 1 : 1 : degree + 1
    
    rightControlPoints(i,:) = Pk_array{i}(end, :);
    
end

rightControlPoints(i,:) = Pk_array{i}(end,:);

rightControlPoints = sortrows(rightControlPoints,1);

end