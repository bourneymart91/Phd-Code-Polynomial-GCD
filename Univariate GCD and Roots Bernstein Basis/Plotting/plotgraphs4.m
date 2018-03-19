
function [] = plotgraphs4(alphas,thetas,residual)

% This function plots graphs of the values of alpha, theta and the
% residual.    



% Graph 1: The graph of the residuals.
figure_name = sprintf([mfilename ' : Residuals' ]);
figure('name',figure_name)

plot(1:1:length(residual),log10(residual),'blue-o','LineWidth',1,...
    'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10);
axis([1,length(residual),-inf,inf])
xlabel('iteration','FontSize',16)
ylabel('log_{10} residual','FontSize',16)
title('Residual','FontSize',14)



%%%%%%%%%%%%%%%%%%%%%%%

% Graph 2: The graph of alpha

% Graph 2: The graph of res_vx
figure_name = sprintf([mfilename ' : Residual v(x)']);
figure('name',figure_name);

        
xalpha=1:1:length(alphas);
plot(xalpha,alphas,'-bs', 'LineWidth',1,'MarkerEdgeColor','b',...
    'MarkerFaceColor','b','MarkerSize',10); 
axis([1,length(alphas),-inf,inf])
xlabel('k','FontSize',16)
ylabel('\alpha','FontSize',16)
title('alpha ','FontSize',14);   
      

        
%%%%%%%%%%%%%%%%%%%%%%%

% Graph 3: The graph of theta.

figure_name = sprintf([mfilename ' : Graphing theta)']);
figure('name',figure_name);


xtheta=1:1:length(thetas);
plot(xtheta,thetas,'-bs', 'LineWidth',1,'MarkerEdgeColor','b',...
    'MarkerFaceColor','b','MarkerSize',10); 
axis([1,length(thetas),-inf,inf])
xlabel('k','FontSize',16)
ylabel('\theta','FontSize',16)
title('theta ','FontSize',14);   



end


    


