set(gca,'FontSize',18)

hLines = findobj(gca,'Type','line')



for k = 1:length(hLines)
    activeLine = hLines(k) 
    
    %set(activeLine,'LineWidth',6)
    %set(activeLine,'markers',1)
    set(activeLine,'Marker','*')
    %set(activeLine, 'MarkerEdgeColor','black')
    %set(activeLine,'MarkerFaceColor','black')
    
end

grid on
box on

myplot = gcf;
filename = myplot.FileName;
[pathstr,name,ext] = fileparts(filename) ;

name1 = strcat(name,'.eps');
saveas(myplot,fullfile(pathstr, name1),'epsc')

name2 = strcat(name,'.fig');
saveas(myplot,fullfile(pathstr, name2),'fig')