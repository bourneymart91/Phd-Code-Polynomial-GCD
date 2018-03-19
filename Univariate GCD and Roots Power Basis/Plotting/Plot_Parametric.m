function [] = Plot_Parametric(at,bt,wt)
% Plot the rational parametric polynomial curve (x(t),y(t))
% where:
%   x(t) = a(t) / w(t)
%   y(t) = b(t) / w(t)

%%
t = -200:0.1:200;


%Note that polyval considers polynomial coefficients highest power first,
%so flip ordering, since my work uses lowest power first.
as = polyval(flipud(at),t);
bs = polyval(flipud(bt),t);
ws = polyval(flipud(wt),t);


figure('name','Parametric Plot')
hold on
plot(as./ws,bs./ws)
hold off

end
