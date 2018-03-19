myplot = gcf;
filename = myplot.FileName;
[pathstr,name,ext] = fileparts(filename) ;


% For Third Size Plots
%label_size = 20;
%tick_size = 1;



% set(gca,'FontSize',tick_size);
% 
% xlbl = get(gca,'Xlabel');
% set(xlbl,{'FontSize'},{label_size});
% 
% ylbl = get(gca,'Ylabel');
% set(ylbl,{'FontSize'},{label_size});
% 
% try
%     zlbl = get(gca,'Zlabel');
%     set(zlbl,{'FontSize'},{label_size});
% catch
%     
% end

name1 = strcat(name,'.eps');
saveas(myplot,fullfile(pathstr, name1),'epsc')

name2 = strcat(name,'.fig');
saveas(myplot,fullfile(pathstr, name2),'fig')



