function [] = MySave(title)


bool_save = 'n';
switch bool_save
    case 'y'
name = sprintf('Outputs/%s',title);
saveas(gcf,name,'epsc')
saveas(gcf,name,'jpeg')
saveas(gcf,name,'png')
saveas(gcf,name,'fig')
    case 'n'
end
end