

function [new_left_pk, new_right_pk] = BezierSubdivide(Pk, degree, c)
% Given a set of control points, subdivide such that two new sets of
% control points are obtained.
%
% Inputs.
%
%
% Pk : Set of original control points.
%
% degree : degree of control points
%
% c : point of subdivision.
%
%
% Outputs.
%
%
% new_left_pk : - The set of control points on the left of c.
%
% new_right_pk : - The set of control points on the right of c.
%


% obtain each set of control points
Pk_array{1} = Pk;

for i = 2 : 1 : degree + 1
    Pk_array{i} = deCasteljau(c,Pk_array{i-1})  ;
end

% From the generated sets of control points obtain set of control points
% for sub-curve P[a,c]
new_left_pk = zeros(degree+1,2);
for i = 1:1:degree+1
    new_left_pk(i,:) = Pk_array{i}(1,:);
end

% obtain set of control points for sub-curve P[c,b]
new_right_pk = zeros(degree+1,2);
for i=1:1:degree+1
    new_right_pk(i,:) = Pk_array{i}(end,:);
end

new_right_pk(i,:) = Pk_array{i}(end,:);

new_right_pk = sortrows(new_right_pk,1);

end