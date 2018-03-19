function PlotAlphas(vAlpha)
%
% Used in SNTLN
% 
% vAlpha : (Vector) 

figure_name = sprintf([mfilename ' : Alpha variation over SNTLN']);
figure('name',figure_name)
hold on
plot(log10(vAlpha),'-s','DisplayName','\alpha')
xlabel('Iteration');
ylabel('log_{10}');
hold off

end