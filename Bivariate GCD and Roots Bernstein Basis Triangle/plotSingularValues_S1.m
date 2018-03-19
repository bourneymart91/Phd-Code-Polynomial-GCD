function [] = plotSingularValues_S1(vSingularValues)



figure()

hold on 

plot(log10(vSingularValues),'-s');

hold off

end