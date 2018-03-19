function t = GetGCDDegree_OneSubresultant(vSingularValues)
% Given the vector of values from either minimum singular values or max:min
% R diagonals.
% Get the rank, where only one subresultant exists.


global SETTINGS

[St,~] = dbstack();
calling_function = St(2).name;


% Only one subresultant
fprintf([mfilename ' : ' calling_function ' : ' 'Only one subresultant exists. \n'])
if(SETTINGS.PLOT_GRAPHS)
    
        figure_name = sprintf([calling_function ' : Singular values of S_{1}']);
        figure('name',figure_name)
        hold on
        title('Singular values of S_{1}')
        plot(log10(vSingularValues))
        hold off
    
end

[deltaSingularValues,~] = Analysis(vSingularValues);

% If the change is smaller than the predefined threshold value, then plot
% is considered 'flat'.

fprintf([mfilename ' : ' calling_function ' : ' sprintf('Threshold : %2.4f \n', SETTINGS.THRESHOLD)]);
fprintf([mfilename ' : ' calling_function ' : ' sprintf('Max Change is Singular Values : %2.4f \n' , deltaSingularValues)])

if deltaSingularValues < SETTINGS.THRESHOLD
    
    % The subresultant is of full rank, in which case t = 0
    t = 0;
    fprintf([mfilename ' : ' calling_function ' : ' 'The only Subresultant S_{1} appears to be full rank. \n']);
    return
    
else % val > threshold
    
    % The subresultant S_{1} is rank deficient, in which case t = 1
    t = 1;
    fprintf([mfilename ' : ' calling_function ' : ' 'The only Subresultant S_{1} appears to be rank defficient \n']);
    return
    
end


end