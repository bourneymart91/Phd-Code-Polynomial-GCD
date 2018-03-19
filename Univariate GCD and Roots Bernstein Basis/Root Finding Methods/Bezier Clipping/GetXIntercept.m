function [x_intercept] = GetXIntercept(convex_hull_vertices)
% Get the intersections between the convex hull and the x axis, return the
% first intercept.
%
% Inputs.
%
%
% convex_hull_vertices : Set of [x,y] paris of vertices of the convex hull.
%                        where the first vertex is included twice.
%
% Outputs.
%
% x_intercept : the first point at which the convex hull crosses the x
%               axis.
%



% Get the number of vertices in the convex hull
[nVerticesCH,~] = size(convex_hull_vertices);

% Note that convex_hull_vertices contains the first vertex twice
nVerticesCH = nVerticesCH - 1;

% Get the set of edges which make up the convex hull
nEdges = nVerticesCH;

% Initialise a vector to store xaxis intercepts
v_x_intercept = [];

for i =1:1:nEdges
    x0 = convex_hull_vertices(i,1);
    y0 = convex_hull_vertices(i,2);
    x1 = convex_hull_vertices(i+1,1);
    y1 = convex_hull_vertices(i+1,2);
    
    if (y0*y1) > 0 % No change of sign
        % Do Nothing
    else % Change of sign
        
        
        % Get x intercept
        m = (y1-y0)./ (x1-x0);
        
        % Get location of intercept on x axis
        location = x1 - (y1./m);
        
        v_x_intercept = ...
            [
            v_x_intercept ;
            location
            ];
        
        
    end
end

if size(v_x_intercept,1) < 1
    
    x_intercept = -1000;
    return;
end

% Get all x intercepts in order
v_x_intercept = sort(v_x_intercept,'ascend');

% Get the first intercept
x_intercept = v_x_intercept(1);
end