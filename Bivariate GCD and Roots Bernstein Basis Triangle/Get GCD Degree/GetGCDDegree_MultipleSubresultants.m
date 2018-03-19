function [t] = GetGCDDegree_MultipleSubresultants(vMetric, limits_k, limits_t,  rank_range)
%
% % Inputs
%
% vMetric : (Vector) 
%
% limits_k : (Int Int)
%
% rank_range : [Float Float]
%
% % Outputs
%
% t : (Int) Degree of the GCD


% Get the name of the function which called this function
[St, ~] = dbstack();
calling_function = St(2).name;

% Get upper and lower bounds
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);


lowerLimit_t = limits_t(1);
upperLimit_t = limits_t(2);


rank_range_low = rank_range(1);
rank_range_high = rank_range(2);

% Get previous change in metric from previous GCD computation
previousDelta = abs(diff(rank_range));


vMetric = sanitize(vMetric);

vDeltaMetric = abs(diff(vMetric));

[maxDelta, indexMaxDelta] = max(vDeltaMetric);

vec_k = lowerLimit_k : 1 : upperLimit_k - 1;

t_pressumed = vec_k(indexMaxDelta);

while ( t_pressumed < lowerLimit_t )

    % Remove
    vDeltaMetric(indexMaxDelta) = []; 
    vec_k(indexMaxDelta) = [];

    % Get new max delta
    [maxDelta, indexMaxDelta] = max(vDeltaMetric);
    
    % Set presumed 
    t_pressumed = vec_k(indexMaxDelta);

    
end

if isempty(t_pressumed)
    t_pressumed = upperLimit_k;
end

% check if the maximum change is significant
fprintf([mfilename ' : ' sprintf('Max change : %2.4f \n', maxDelta)]);
fprintf([mfilename ' : ' sprintf('Previous Delta : %2.4f \n', previousDelta )]);



if abs(maxDelta) < 2 %(0.55 * previousDelta)
    
    
    fprintf([calling_function ' : ' mfilename ' : ' 'Polynomials either coprime or GCD = g(x) \n' ])
    
    
    % Change in Singular values is not significant so check if all
    % subresultants are rank deficient or full rank
    
    % Get the average of the minimum singular values
    avgMetricValue = mean(vMetric);
    
    fprintf([calling_function ' : ' mfilename ' : ' sprintf('Average Singular Value : %e \n',avgMetricValue) ])
    
    if avgMetricValue > mean(rank_range)
       % All Minimum singular values are below threshold so, all 
       % subresultants are rank deficient. deg(GCD) = 0
       fprintf([calling_function ' : ' mfilename ...
           ' : ' ...
           'All Subresultant matrices are rank deficient \n' ])
       t = upperLimit_k;
       
    else 
        % All minimum singular values are above threshold so all
        % subresultants are full rank. deg(GCD) = min(m,n)
       fprintf([calling_function ' : ' mfilename ' : ' 'All Subresultants are rank deficient : GCD = g(x) \n' ])
       t = 0;
       
    end

else
    % change is significant
    t = t_pressumed;

    
end

end