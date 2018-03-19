
function [roots, root_count, reached_depth] = ...
    FindRoots(CP, degree, curr_depth, W_DEGREE, root_count, reached_depth)
%
% Inputs
%
% CP : (Matrix) Each row contains the coordinate pair of the control points
% [x_{i},y_{i}] of f(x)
%
% degree : (Int) Degree of the polynomial f(x)
%
% curr_depth - current depth of iteration
%
% W_Degree - Degree of control points
%
% root count - Number of roots obtained to this depth.
%
% reached depth - maximum reached depth
%
% t -
%

printSpacing(curr_depth)

fprintf('Inspecting interval %f - %f \n', CP(1,1), CP(end,1))

% If the current depth is greater than reached depth so far, increment by 1
if curr_depth > reached_depth
    
    reached_depth = curr_depth;
    
end

global MAXDEPTH

% Initialise empty set of roots
roots = [];

% Get number of crossings of control polygon and x axis.
nCrossings = CrossingCount(CP);

% get the number of crossings of the control polygon.
switch(nCrossings)
    case 0 % No solutions in the interval
        
        printSpacing(curr_depth)
        fprintf('No Crossings in the interval %f - %f \n', CP(1,1), CP(end,1))
        
    case 1
        % If number of crossings is 1, then there is a unique solution in
        % the interval
        
        printSpacing(curr_depth)
        fprintf('One Crossing in the interval %f - %f \n', CP(1,1), CP(end,1))
        
        if (curr_depth >= MAXDEPTH) % Current depth is equal to or greater than the maximum depth
            
            
            
            % Calculate an approximate intercept, in the middle of the
            % interval.
            t_new = (CP(1,1) + CP(W_DEGREE+1,1)) ./2;
            
            printSpacing(curr_depth)
            fprintf('Maximum Depth Reached \n')
            printSpacing(curr_depth)
            fprintf('** Root Approximated at %f \n',t_new)
            
            % Add the approximate root to the set of roots
            roots = [roots;t_new];
            
            % Increment the root count.
            root_count = root_count + length(t_new);
            
            
            
        elseif (ControlPolygonFlatEnough(CP, degree))
            % if the control polygon is considered flat enough, that the
            % curve which it bounds may be considered a line. Calculate
            % the x intercept
            
            printSpacing(curr_depth)
            fprintf('Control Polygon Flat Enough\n')
                        
            % Calculate the new root
            t_new = ComputeXIntercept(CP, degree);
            
            printSpacing(curr_depth)
            fprintf('** Root at %f \n', t_new)
            
            
            % Add the root to the list of roots.
            roots = [roots; t_new];
            
            % Increment the root count
            root_count = root_count + length(t_new);
            
            
        else
            
            printSpacing(curr_depth)
            fprintf('Control Polygon not flat enough for interval %f - %f \n', CP(1,1), CP(end,1))
            
            printSpacing(curr_depth)
            fprintf('Subdivide Interval \n')
            
            % If the control polygon isn't flat enough. subdivide the
            % control polygon and look at the left and right control
            % polygons individually.
            
            % Get start point
            a = CP(1,1);
            
            % Get end point
            b = CP(degree+1,1);
            
            % Get midpoint
            c = a + ((b-a) /2);
            
            % Get a set of control points for the left and right halves of
            % the interval.
            [leftControlPoints, rightControlPoints] = BezierSubdivide(CP, degree, c);
            
            
            % Perform find roots on the left partition, while incrementing
            % the current depth by one.
            [t_left,~, reached_depth]  = FindRoots(leftControlPoints, degree, curr_depth+1, W_DEGREE, root_count, reached_depth);
            
            % Add list of roots returned by findRoots() on the left
            % partition, to the list of all roots.
            roots = [roots; t_left];
            
            % Increment the root count.
            root_count = root_count + length(t_left);
            
            
            % Perform findRoots() on the right set of control points
            [t_right,~,reached_depth] = FindRoots(rightControlPoints, degree, curr_depth+1, W_DEGREE, root_count, reached_depth);
            
            % Add list of roots returned by findRoots() on the right
            % partition, to the list of all roots.
            roots = [roots; t_right];
            
            % Increment root count.
            root_count = root_count + length(t_right);
            
        end
        
    otherwise
        % If the number of crossings is not 0 or 1, then number of crossings is greater than 1.
        % The interval contains more than one root.
        
        
        
        if (curr_depth >= MAXDEPTH) % Current depth is equal to or greater than the maximum depth
            
            printSpacing(curr_depth);
            fprintf('Maximum Depth Reached \n')
            
            % Calculate an approximate intercept, in the middle of the
            % interval.
            t_new = (CP(1,1) + CP(W_DEGREE+1,1)) ./2;
            
            % Add the approximate root to the set of roots
            roots = [roots; t_new];
            
            % Increment the root count.
            root_count = root_count + length(t_new);
        else
            
            % Split the interval
            printSpacing(curr_depth)
            fprintf('Control polygon has multiple crossings so subdivide\n')
            
            printSpacing(curr_depth)
            fprintf('Subdivide Interval \n')
            
            % Get Start Point
            a = CP(1,1);
            
            % Get End Point
            b = CP(degree+1,1);
            
            % Get Midpoint
            c = a + ((b-a) /2);
            
            % Get a set of control points for the left and right halves of
            % the interval.
            [leftControlPoints, rightControlPoints] = BezierSubdivide(CP, degree, c);
            
            
            % Perform find roots on the left partition, while incrementing
            % the current depth by one.
            [t_left,~,reached_depth]  = FindRoots(leftControlPoints, degree, curr_depth+1, W_DEGREE, root_count, reached_depth);
            
            % Add list of roots returned by findRoots() on the left
            % partition, to the list of all roots.
            roots = [roots; t_left];
            
            % Increment root count.
            root_count = root_count + length(t_left);
            
            % Perform findRoots() on the right set of control points
            [t_right,~,reached_depth] = FindRoots(rightControlPoints, degree, curr_depth+1, W_DEGREE, root_count, reached_depth);
            
            % Add any calculated roots to the list of all roots.
            roots = [roots; t_right];
            
            % Increment root count.
            root_count = root_count + length(t_right);
        end
end


end


function printSpacing(curr_depth)
for i = 1:1:curr_depth
    fprintf('\t')
end
end