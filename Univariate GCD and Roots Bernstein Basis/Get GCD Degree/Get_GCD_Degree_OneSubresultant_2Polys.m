function t = Get_GCD_Degree_OneSubresultant_2Polys(vMetric, rank_range)
% Given the vector of values from either minimum singular values or max:min
% R diagonals.
% Get the rank, where only one subresultant exists.
%
% Inputs
%
% vSingularValues : (Vector) Vector of rank revealing metric
%
% rank_range : [Float Float] 
%
% Outputs.
%
% t : (Int) Computed degree of the GCD


global SETTINGS


% 
rank_range_low = rank_range(1);
rank_range_high = rank_range(2);
previousDelta = abs(diff(rank_range));

% Get the calling function
[St,~] = dbstack();
calling_function = St(2).name;

% Only one subresultant
fprintf([calling_function ' : ' mfilename ' : ' 'Only one subresultant exists. \n'])

% Plot the singular values
if(SETTINGS.PLOT_GRAPHS)
    figure_name = sprintf([calling_function ' : Singular values of S_{1}']);
    figure('name',figure_name)
    hold on
    title('Singular values of S_{1}')
    plot(vMetric)
    hline(rank_range_low,'r')
    hline(rank_range_high,'r')
    hline(mean(rank_range),'b')
    
    hold off
    
end

average_singularValues = mean(vMetric);

[maxDelta,~] = Analysis(vMetric);

% If the change is smaller than the predefined threshold value, then plot
% is considered 'flat'.
if maxDelta < (0.9 * previousDelta)
    
    
    % All full rank or all rank deficient
    
    
    if average_singularValues > rank_range_high
        
        % all full rank
        t = 0;
    elseif average_singularValues < rank_range_low
            % all rank deficient
            t = 1;
            
    else
        % Average is somewhere in the middle
        t = 1;
    end
    
    fprintf([mfilename ' : ' calling_function ' : ' sprintf('t = %i \n',t)]);
    return
    
else % val > threshold
    
    % The subresultant S_{1} is rank deficient, in which case t = 1
    t = 1;
    fprintf([mfilename ' : ' calling_function ' : ' 'The only Subresultant S_{1} appears to be rank deficient \n']);
    fprintf([mfilename ' : ' calling_function ' : ' sprintf('t = %i \n',t)]);
    return
    
end
end