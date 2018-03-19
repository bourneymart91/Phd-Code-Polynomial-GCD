
function [Pk_1] = deCasteljau(c, Pk_0)
% Obtain the next set of control points, given a set of control points Pk_0
% in the deCasteljau subdivision process.
% Uses deCasteljau algorithm found in 'Applied Geometry for computer
% graphics and CAD' - Page 151
%
% % Inputs
%
% c : (Float) subdivision point on horizontal axis
%
% Pk_0 : (Matrix) Each row contains the coordinates of a control point

% get first control poitn
a = Pk_0(1,1);

% get last control point
b = Pk_0(end,1);


T = (c-a)./(b-a);


% Build set of control points
Pk_1 = [];
% for
for i = 1:1:length(Pk_0)-1
    %Pk_1(i)
    aa =  (1-T).*Pk_0(i,:) + T.*Pk_0(i+1,:) ;
    Pk_1 = [Pk_1 ; aa];
    
end



end