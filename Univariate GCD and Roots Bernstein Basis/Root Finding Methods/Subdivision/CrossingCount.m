function nCrossings = CrossingCount(CP)
% Get number of crossings between a given control polygon and the X axis.
%
% % Inputs
%
% CP : (Matrix) Each row contains a coordinate pair of the control 
% points [x,y]
%
% Outputs.
%
% nCrossings : (Int) Number of crossings of the control polygon and the x
% axis



% Initialise old sign to be the value of y for the first control point.
old_sign = CP(1,2);

% Initialise number of crossings = 0.
nCrossings = 0;

% For each control point get its sign, and the sign of the previous control
% point, and compare them.
for i = 1 : 1 : size(CP, 1) - 1
    
    % Assign the new sign to be y value for current control point
    new_sign = CP(i + 1, 2);
    
    % If there is a change of sign y0*y1 < 0, so incrememnt number of
    % crossings.
    if (new_sign * old_sign) < 0
        
        % if change of sign, then increment number of crossings
        nCrossings = nCrossings + 1;
        
    end
    
    % Set old sign to be current iterations new sign, then start next
    % iteration.
    old_sign = new_sign;
end

end
