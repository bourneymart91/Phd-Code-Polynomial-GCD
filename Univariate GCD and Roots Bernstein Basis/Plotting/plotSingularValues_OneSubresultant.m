function [] = plotSingularValues_OneSubresultant(vSingularValues, k)


figure_name = sprintf('Singular values of S%i',k);
figure('name',figure_name)

hold on

x_vec = 1:1:length(vSingularValues);

plot(x_vec, log10(vSingularValues),'-s')

xlabel('$i$','Interpreter','latex')
ylabel('$\log_{10}\left( \sigma_{k,i} \right)$','Interpreter','latex');

hold off

% Display of graphs
myplot = gca;

myval_side = 0.12;
myval_base = 0.10;

%set(gca,'FontSize',20)
%set(gca,'FontName','Helvetica')

set(myplot, 'Position', [ myval_side myval_base 0.98 - myval_side 0.98 - myval_base])

grid on
box on

end