function plotRowNorm_degreeTotal(arr_R1_RowNorm,  limits_k, limits_t, rank_range)
%
% % Inputs
%
% arr_R1_RowNorm : (Array of Vectors)
%
% limits_k : My Limits on degree of GCD computation
%
% limits_t : Actual limits on degree of GCD computation

% Get my upper and lower limit
lowerLimit_k = limits_k(1);
upperLimit_k = limits_k(2);


% Get number of subresultant matrices
nSubresultants = upperLimit_k - lowerLimit_k + 1;


figure_name = sprintf([mfilename ' : ' 'Row Norms']);
figure('name',figure_name);
hold on

for i = 1:1:nSubresultants
    
    vR1RowNorm = arr_R1_RowNorm{i};
    
    k = lowerLimit_k + (i-1);
    
    % Produce a vector of i
    vec_k = k.* ones(length(vR1RowNorm),1);
    
    plot(vec_k, log10(vR1RowNorm),'*');
    
end

vline(limits_t, {'-r','-r'});
hline([rank_range mean(rank_range)],{'-r','-r','-r'});


hold off


end