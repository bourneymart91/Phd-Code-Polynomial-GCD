function PlotThetas(vTheta)
%
% % Inputs
%
% vTheta : vector of theta values for each iteration

figure_name = sprintf([mfilename ' : Theta variation over SNTLN']);
figure('name',figure_name)
hold on
plot(log10(vTheta),'-s','DisplayName','\theta')
xlabel('Iteration');
ylabel('log_{10}');
hold off

end