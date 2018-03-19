function [t] = GetGCDDegree_MultipleSubresultants(vMetric, limits_k, limits_t, rank_range)
% Get the problem type, dependent on the vector of singular values from the
% series s_{k}
%
% % Inputs
%
% vMinimumSingularValues : (Vector)
%
% limits_k : [Int Int] 
%
% limits_t : [Int Int]
%
% rank_range : [Float Float]
%
%
% Get the type of problem.
% Problem Type.
% Singular      : All Subresultants S_{k} are Singular, and rank deficient
% NonSingular   : All Subresultants S_{k} are Non-Singular, and full rank
% Mixed         : Some Subresultants are Singular, others are Non-Singular.

% Get the name of the function which called this function
[St,~] = dbstack();
calling_function = St(2).name;

% Set upper and lower limits of the degree of the GCD.
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);



% Set global variables
global SETTINGS

% Get the index of the largest change in minimum singular values
[max_delta, indexMaxChange] = Analysis(vMetric);

MessageToConsole( sprintf('Largest Change in log of Singular Values : %4.5e' ,abs(max_delta)));
MessageToConsole( sprintf('Threshold : %4.5e', SETTINGS.THRESHOLD));

% Check if the delta is significant
if (abs(max_delta) < SETTINGS.THRESHOLD)
    
    bool_significant_change = false;
    
else
    
    bool_significant_change = true;
    
end


% if the largest of the changes is below a threshold, then the largest
% change is not significant, and cannot determine the degree of the GCD. so
% all subresultants are either singular or nonsingular. The GCD is either 0
% or k where k is the upper limit.
if  (bool_significant_change == false)
    

    % % Rank Deficient = Singular = GCD = maximum
    % % Full Rank = Non-Singular => GCD = 0
    % maxChange is insignificant
    % Get the average minimum singular value
    avgMinSingularValue = (mean(vMetric));
    
    MessageToConsole( sprintf('Threshold to dermine rank : %e \n', SETTINGS.THRESHOLD_RANK));
    MessageToConsole( sprintf('Average of Singular values: %2.4f \n',avgMinSingularValue) );
    
  
    
    if  avgMinSingularValue < mean(rank_range)
        % If all singular values are close to zero, then the Sylvester matrices
        % are all rank deficient, and are all singular
        % gcd is min(m,n)
        t = upperLimit_k;
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('All Rank Deficient \n')]);
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('t = %i \n',t)]);
    
        
    elseif lowerLimit_k > 1
            t = lowerLimit_k;
            fprintf([mfilename ' : ' calling_function ' : ' sprintf('t = %i \n',t)]);
        
    else
        t = 0;
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('All Full Rank \n')]);
        fprintf([mfilename ' : ' calling_function ' : ' sprintf('t = %i \n',t)]);
        
    end
else
    
    t = lowerLimit_k + indexMaxChange - 1;
    % maxChange is signifcant
    fprintf([mfilename ' : ' calling_function ' : ' 'min < Deg(GCD) < max \n'])
    fprintf([mfilename ' : ' calling_function ' : ' sprintf('t = %i \n',t)])
    
end



end

function [] = MessageToConsole(str)
[St,~] = dbstack();
calling_function = St(3).name;


fprintf([mfilename ' : ' calling_function ' : ' str '\n'])
end