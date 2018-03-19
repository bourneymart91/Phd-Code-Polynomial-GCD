


function PlotResiduals(vCondition)


figure_name = sprintf([mfilename ' : ' 'Residuals']);
figure('name',figure_name)
hold on
xlim([1 +inf]);
plot(log10(vCondition),'-s');
hold off
        
end

