function [t] = GetGCDDegree_OneSubresultant(Sk)
% % Inputs
%
% Sk : (Matrix) kth Sylvester subresultant matrix
%
% % Outputs
%
% t : (Int) Degree of GCD

global SETTINGS

fprintf('Only one subresultant exists. Check if near rank deficient \n')

vSingularValues = svd(Sk);


% Plot all singular values of S_{1}(f,g)
figure('name','GetDegreeTotal - SVD')
hold on
plot(log10(vSingularValues),'-s')
hold off


% Get deltas
vDelta_MinSingularValues = abs(diff(vSingularValues));

% Get maximum delta
max_change = max(log10(vDelta_MinSingularValues));

if max_change < SETTINGS.THRESHOLD
    %
    fprintf([mfilename 'Change in Singular values is not significant \n']);
    t = 0;
    
    
else
    fprintf([mfilename 'Change in Singular values is signficant \n']);
    t = 1;
end

fprintf('Degree of GCD : %i \n',t)



end