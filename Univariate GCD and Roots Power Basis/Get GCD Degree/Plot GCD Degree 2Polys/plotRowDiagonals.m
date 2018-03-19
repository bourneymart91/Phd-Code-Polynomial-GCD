
function plotRowDiagonals(arr_RowDiagonals, myLimits, limits)

% % Inputs
%
% arr_RowDiagonals : Array containing the diagonals of the matrices R_{k}, 
% from the QR decomposition of S_{k} for k = lower_lim:upper_lim 
%
% myLimits : [(Int) (Int)] : Limits on the number of subresultants computed 
%
% k_limits : [(Int) (Int)] : Limits on the possible values of k
%
%

% Get my Limits
myLowerLimit = myLimits(1);
myUpperLimit = myLimits(2);

% Get number of subresultants
nSubresultants = myUpperLimit - myLowerLimit + 1;

% Get limits on degree of GCD
lowerLimit = limits(1);
upperLimit = limits(2);

% Plot graph of norms of each row (N) from the qr decompostion of each S_{k}
figure_name = sprintf('%s : Row Sum Norm',mfilename);
figure('name',figure_name)
hold on

for i = 1 : 1 : nSubresultants
    
    k = lowerLimit + (i-1);

    % Get set of diagonals
    vec_RowDiags = arr_RowDiagonals{i};
    vec_k = k.*ones(length(vec_RowDiags));
    
    plot(vec_k, log10(vec_RowDiags) ,'*')
end

xlabel('k')

ylabel(sprintf('Diagonals of R1 from QR decomposition of S_{k}'))
title(sprintf('Diagonals of R1 from QR decomposition of S_{k}'));

xlim([1 upperLimit]);
vline(lowerLimit,'b','');
vline(upperLimit,'b','');

hold off
end
