function [t] = GetGCDDegree_MultipleSubresultants(vMetric, limits_k, limits_t, rank_range)
%
% Given a metric - provided by one of the GetGCDDegree functions, determine
% the degree of the GCD.
%
% % Inputs
%
% vMetric : vector of values from the metric used to determine whether each
% sylvester subresultant matrix is full rank or rank deficient.
%
% limits_k : [Int Int] index of first and last
% sylvester subresultant whose rank revealing metric is computed.
% 
% limits_t : [Int Int] 
%
% rank_range : [Float Float]
%

% Get the name of the function which called this function
[St,~] = dbstack();
calling_function = St(2).name;

% Set Global Variables.

% Get upper and lower bound on degree of GCD
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(2);

% Sanitize the metric.
vMetric = Sanitize(vMetric);


% Get the differences between minimum singular value S_{k} and S_{k+1}
vDeltaMetric = abs(diff(vMetric));
v_k = lowerLimit_k : 1 : upperLimit_k-1;

% Get the maximum change (on log scale) in minimum singular values.
[maxDelta, index] = max((vDeltaMetric));
k = v_k(index);


if ( k < lowerLimit_t)
   
    % Remove the max value from the matrix
    vDeltaMetric(index) = [];
    v_k(index) = [];
    
    [maxDelta, index] = max(vDeltaMetric);
    
    k = v_k(index);
    
end


% Get previous delta
previousDelta = abs(diff(rank_range));

fprintf([mfilename ' : ' sprintf('Delta : %4.5e \n' ,abs(maxDelta))]);
fprintf([mfilename ' : ' sprintf('Delta : %4.5e \n' ,abs(previousDelta))])



% check if the maximum change is significant
if abs(maxDelta) < 0.5*previousDelta
    
    
    if mean(vMetric) < mean(rank_range)
        % All singular values are sufficiently small for the subresultants
        % to all be rank deficient
        t = upperLimit_k;
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('All Subresultant matrices are singular (rank deficient) \n')]);
        
    else
        % All singular values are large, and indicate all subresultants are
        % full rank
        t = 0;
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('All subresultant matrices are non-singular (full rank) \n')]);
        
    end
    
    fprintf([mfilename ' : ' calling_function ' : ' sprintf('Degree of the GCD : %i \n',t)]);
    
    
else
    % max_change is sufficiently large to indicate degree of GCD.
    t = lowerLimit_k + index - 1;
    fprintf([mfilename ' : ' calling_function ' : ' sprintf('Degree of the GCD : %i \n',t)])
    
end