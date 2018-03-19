function [maxDelta, indexMaxDelta] = Analysis(vMetric)
%
% % Inputs
%
% vMetric : (Vector) Vector of rank revealing metric
%
% % Outputs
%
% maxDelta : Maximum Difference in the vector of metric values
%
% indexMaxDelta : Index of maximum change

    vMetric = sanitize(vMetric);
    % % Analyse Max:Min Row Norms for each subresultant
    
    % Get the change in the ratios from one subresultant to the next.
    vDeltaMetric = abs(diff((vMetric)));
    
    % Get the maximum change in rowsum ratio and its index
    [maxDelta,indexMaxDelta] = max(vDeltaMetric);
    
end