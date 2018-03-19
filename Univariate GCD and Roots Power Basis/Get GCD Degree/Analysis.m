
function [max_Delta_MaxMin_Diag_R, indexMaxChange_Ratio_DiagsR] = Analysis(vMetric)

% Get Degree by max:min Diagonals

vMetric = sanitize(vMetric);

% Get the change in the ratios of diagonal elements from one subresultant
% to the next.
vDelta_MaxMin_Diag_R = abs(diff(vMetric));

% Get the maximum change in diag ratio and its index
[max_Delta_MaxMin_Diag_R, indexMaxChange_Ratio_DiagsR] = max(vDelta_MaxMin_Diag_R);


end