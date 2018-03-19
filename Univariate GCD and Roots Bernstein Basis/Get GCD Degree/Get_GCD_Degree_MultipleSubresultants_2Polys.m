function [t] = Get_GCD_Degree_MultipleSubresultants_2Polys(vMetric, limits_k, limits_t, rank_range)
% Get the problem type, dependent on the vector of singular values from the
% series s_{k}
%
% Get the type of problem.
% Problem Type.
%   Singular      : All Subresultants S_{k} are Singular, and rank deficient
%   NonSingular   : All Subresultants S_{k} are Non-Singular, and full rank
%   Mixed         : Some Subresultants are Singular, others are Non-Singular.
%
% % Inputs
%
% vMetric : (Vector) vector of values used to determine degree of GCD. These values
% may be 
%       'Minimum Singular values of S_{k}'
%       'Min/Max Row Diagonals of R_{k}'
%       'Min/Max Row Norms of R_{k}'
%
% limits_k : [(Int) (Int)] Range of k values for which the Sylvester subresultant matrix
%   S_{k} is constructed
%
% limits_t : [Int Int] 
%
% rank_range : [(Float) (Float)] 



% Get the function which called this function.
[St,~] = dbstack();
calling_function = St(2).name;

lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);

lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(2);

rank_range_low = rank_range(1);
rank_range_high = rank_range(2);

previousDelta = abs(diff(rank_range));

% Get the maximum change in singular values and the index at which the
% maximum change occured.
vMetric = sanitize(vMetric);

% Get the change in the ratios from one subresultant to the next.
vDeltaMetric = abs(diff((vMetric)));

% Get the maximum change in rowsum ratio and its index
[maxDelta, indexMaxDelta] = max(vDeltaMetric);

% Get the vector of all possible values of t
vec_k = lowerLimit_k : 1 : upperLimit_k - 1;

% Get the presumed degree of the GCD (based only on maximum delta)
t_pressumed =  vec_k(indexMaxDelta);

% If presumed degree is less than the lower limit, remove it from the
% vector of possible degrees, and remove corresponding delta from the
% vector or deltas. Look for the next smallest delta, until the index of
% the maximum delta is within the upper and lower bounds of t.

while ( t_pressumed < lowerLimit_t )

    % Remove
    vDeltaMetric(indexMaxDelta) = []; 
    vec_k(indexMaxDelta) = [];

    % Get new max delta
    [maxDelta, indexMaxDelta] = max(vDeltaMetric);
    
    % Set presumed 
    t_pressumed = vec_k(indexMaxDelta);

    
end
    


display([mfilename ' : ' sprintf('Previous Delta : %2.4f', previousDelta)]);
display([mfilename ' : ' sprintf('Current Delta : %2.4f', maxDelta)]);

% Two condition : 
% 1. Is the maxDelta large enough
% 2. Is the index of the maxDelta within the upper and lower bounds of t

% Check conditon 1.
if  abs(maxDelta) < (0.5 * previousDelta)
    
    fprintf('Delta is insignificant \n');
    avgMetricValue = mean(vMetric);
    
    
    display(avgMetricValue)
    
    
    if  avgMetricValue < mean(rank_range)
        
        % If all singular values are close to zero, then rank deficient, degree of
        % gcd is min(m,n)
        t = upperLimit_k;
        
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('All Subresultants are Rank Deficient \n')])
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('t = %i \n',t)])
        
    else
        % If all Rank metric values are not close to zero, then the matrices are all of full rank, degree
        % of gcd is 0
        
        t = 0;
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('All subresultants Full Rank \n')])
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('t = %i \n',t)])
        
    end
    
else
    
    % To do - Add code to make sure that the max change is within the
    % bounds
    
    t = t_pressumed;
        
    
    
    fprintf([mfilename ' : ' calling_function ' : ' 'Mixed \n'])
    fprintf([mfilename ' : ' calling_function ' : ' sprintf('t = %i \n',t)])
    
    
end

end