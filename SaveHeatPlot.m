ax = gca;

%ax.ColorbarVisible = 'on'

%ax.CellLabelFormat = '%i'


ax.OuterPosition = [0 0 1 1];
myplot = gcf;
filename = myplot.FileName;
[pathstr,name,ext] = fileparts(filename) ;

name1 = strcat(name,'.eps');
saveas(myplot,fullfile(pathstr, name1),'epsc')

name2 = strcat(name,'.fig');
saveas(myplot,fullfile(pathstr, name2),'fig')