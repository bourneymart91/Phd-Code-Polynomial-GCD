function SubPlot
% -----------------------------------------------------------
% General usage: This function inserts various MATLAB figure (.fig) files
%                into one figure with multiple subplots
%
% Notes:
%      All desired .fig files should contain only one figure inside it (no
%      subplot inside the .fig files).
%      All desired .fig files should be defined in the 2D space.
%
% Author: Farhad Sedaghati, PhD student at The University of Memphis
% Contact: <farhad_seda@yahoo.com>
% Written: 04/08/2015
% Updated: 06/21/2015
% Updated: 07/14/2015
% Updated: 08/13/2015
%
% -----------------------------------------------------------





% get handle to axes of figure
ax = findobj(gcf,'type','axes');

nFigures = length(ax);


K=input('How many figures in a row do you want to have? \n');
if isempty(K)
    K=2;
end

N = nFigures;

figure;
for i=1:N
    % create and get handle to the subplot axes
    s(i) = subplot(ceil(N/K),K,i); 
    
    % get handle to all the children in the figure
    aux=get(ax(i),'children');
    
    for j=1:size(aux)
        fig(i) = aux(j);
        copyobj(fig(i),s(i)); 
        hold on
    end
    % copy children to new parent axes i.e. the subplot axes
    xlab=get(get(ax(i),'xlabel'),'string');
    ylab=get(get(ax(i),'ylabel'),'string');
    tit=get(get(ax(i),'title'),'string');
    xlabel(xlab);ylabel(ylab);title(tit);
end








