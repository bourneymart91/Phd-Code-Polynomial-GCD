
function plotMinimumResiduals(vMinimumResidual, myLimits, limits)
%
% % Inputs
%
% vMinimumResidual : Vector containing minimum residuals
%
% myLimits :
%
% limits : [lowerLimit : upperLimit]

lowerLimit = limits(1);
upperLimit = limits(2);

myLowerLimit = myLimits(1);
myUpperLimit = myLimits(2);

% % Plot Residuals
figure_name = sprintf('%s : Residuals',mfilename);
figure('name',figure_name);
hold on
x_vec = myLowerLimit: 1: myUpperLimit;
plot(x_vec,log10(vMinimumResidual),'-s','DisplayName','SVD')
ylabel('log r(k)')
xlabel('k')
legend(gca,'show');
vline(lowerLimit);
vline(upperLimit);
hold off
end
