function [] = PlotCoefficients(fx,fw,name)
% Plot the coefficients of the polynomial f(x) and f(w).
%
% Inputs.
% 
% fx : Coefficients of f(x)
%
% fw : Coefficients of f(w)
%
% name : name of function

% //TODO

str = sprintf('Coefficients of %s',name);
figure('name',str)
label1 = sprintf('%s(x)',name);
label2 = sprintf('%s(\\omega)',name);
plot((fx),'-s','DisplayName',label1)
hold on
plot((fw),'-o','DisplayName',label2)
legend(gca,'show');
hold off

end


