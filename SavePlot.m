myplot = gcf;
filename = myplot.FileName;
[pathstr,name,ext] = fileparts(filename) ;

name1 = strcat(name,'.eps');
saveas(myplot,fullfile(pathstr, name1),'epsc')

name2 = strcat(name,'.fig');
saveas(myplot,fullfile(pathstr, name2),'fig')

name3 = strcat(name,'.png');
saveas(myplot,fullfile(pathstr, name3),'png')