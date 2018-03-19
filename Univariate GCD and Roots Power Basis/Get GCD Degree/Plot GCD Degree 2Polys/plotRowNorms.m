
function plotRowNorms(arr_RowNorms, myLimits, limits)
% 
% % Inputs
%
% arr_RowNorms :
%
% myLimits :
%
% limits : 
%
% 


%
myLowerLimit = myLimits(1);
myUpperLimit = myLimits(2);

nSubresultants = myUpperLimit - myLowerLimit + 1;

lowerLimit = limits(1);
upperLimit = limits(2);

figure_name = sprintf('%s : Diag Norm',mfilename);
figure('name',figure_name)
hold on
for i = 1:1:nSubresultants
    
    k = myLowerLimit + (i-1);
    
    % get vector of row norms
    vec_RowNorms = arr_RowNorms{i};
    vec_k = k.*ones(length(vec_RowNorms));
    
    plot(vec_k, log10(vec_RowNorms),'*');

end

hold off
xlabel('k')
ylabel('log10 Row Norm of R1 from QR decomposition of S_{k}')
title(sprintf('log10 Row Norm of R1 from the QR decomposition of each subresultant S_{k}'));
vline(lowerLimit,'b','');
vline(upperLimit,'b','');
hold off
end





