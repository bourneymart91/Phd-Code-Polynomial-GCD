
function [] = plotgraphs3(res_ux,res_vx,res_uw,res_vw)
% This function plots graphs of the residuals of the equations

% (lambda ux)(dx) = fx    and    (mu vx)(dx) = gx

% whose residuals are res_ux and res_vx respectively, and the
% residuals of the equations

% (lambda uw)(dw) = fw    and    (mu vw)(dw) = gw

% whose residuals are res_uw and res_vw respectively,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                           Inputs

% res_ux

% res_vx

% res_uw

% res_vw





% Graph 1
% The graph of res_ux

figure_name = sprintf([mfilename ' : Residual u(x)']);
figure('name',figure_name);
x = 1:length(res_ux);  % the graph of res_ux
plot(x,log10(res_ux),'-rs','LineWidth',1, 'MarkerEdgeColor',...
    'r','MarkerFaceColor','r','MarkerSize',10);          
axis([1,length(res_ux),-inf,inf]);
xlabel('iteration','FontSize',16);
ylabel('log_{10} residual','FontSize',16);
title('Residuals from u(y)d(y)=f(y)','FontSize',12);


%%%%%%%%%%%%%%%%%%%%%%%%%

% Graph 2: The graph of res_vx
figure_name = sprintf([mfilename ' : Residual v(x)']);
figure('name',figure_name);

x = 1:length(res_vx);  % the graph of res_ux
plot(x,log10(res_vx),'-bs','LineWidth',1, 'MarkerEdgeColor',...
    'b','MarkerFaceColor','b','MarkerSize',10);          
axis([1,length(res_vx),-inf,inf]);
xlabel('iteration','FontSize',16);
ylabel('log_{10} residual','FontSize',16);
title('Residuals from v(y)d(y)=g(y)','FontSize',12);


%%%%%%%%%%%%%%%%%%%%%%%%%

% Graph 3: The residuals res_ux and res_vx on the same graph.

% Graph 2: The graph of res_vx
figure_name = sprintf([mfilename ' : Residual u(x) and v(x)']);
figure('name',figure_name);


x = 1:length(res_ux);  % the graph of res_ux
plot(x,log10(res_ux),'-rs','LineWidth',1, 'MarkerEdgeColor',...
    'r','MarkerFaceColor','r','MarkerSize',10);          
axis([1,length(res_ux),-inf,inf]);
xlabel('iteration','FontSize',16);
ylabel('log_{10} residual','FontSize',16);
title('Residuals from u(y)d(y)=f(y)(red) and v(y)d(y)=g(y) (blue)',...
    'FontSize',12);
        
% The graph of res_vw.
hold on 
plot(x,log10(res_vx),'-bo','LineWidth',1, 'MarkerEdgeColor',...
    'b','MarkerFaceColor','b','MarkerSize',10);     
hold off
        

%%%%%%%%%%%%%%%%%%%%%%%%%

% Graph 4: The residuals res_uw and res_vw.

% Graph 2: The graph of res_vx
figure_name = sprintf([mfilename ' : Residual u(\omega) and v(\omega)']);
figure('name',figure_name);


x = 1:length(res_uw);  % the graph of res_ux
plot(x,log10(res_uw),'-rs','LineWidth',1, 'MarkerEdgeColor',...
    'r','MarkerFaceColor','r','MarkerSize',10);          
axis([1,length(res_uw),-inf,inf]);
xlabel('iteration','FontSize',16);
ylabel('log_{10} residual','FontSize',16);
title('Residuals from u(w)d(w)=f(w) (red) and v(w)d(w)=g(w) (blue)',...
    'FontSize',12);
        
% The graph of res_vw.
hold on 
plot(x,log10(res_vw),'-bo','LineWidth',1, 'MarkerEdgeColor',...
    'b','MarkerFaceColor','b','MarkerSize',10);     
hold off
        


end     


    


