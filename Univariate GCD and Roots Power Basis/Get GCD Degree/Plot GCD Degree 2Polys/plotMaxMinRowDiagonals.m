
function plotMaxMinRowDiagonals(vMaxDiagR1,vMinDiagR1, myLimits, limits)
%
% % Inputs
%
% vMaxDiagR1 :
%
% vMinDiagR1 : 
%
% limits :

myLowerLimit = myimits(1);
myUpperLimit = myimits(2);

lowerLimit = limits(1);
upperLimit = limits(2);

x_vec = myLowerLimit : 1 : myUpperLimit;

% Plot Graph of ratio of max : min element of the diagonal elements of R1 from the QR decompositions.
figure_name = sprintf('%s : Max:min Row Diagonals',mfilename);
figure('name', figure_name)

vRatio_MaxMin_Diagonals_R = vMinDiagR1./vMaxDiagR1;
plot(x_vec,log10(vRatio_MaxMin_Diagonals_R),'red-s');
vline(lowerLimit,'b','');
vline(upperLimit,'b','');
hold on
legend('Max:Min diag element of subresultant S_{k}');
title('Max:Min diagonal elements of R1 from the QR decomposition of S_{k} (Original)');
ylabel('log_{10} max:min diag element')
hold off
end
