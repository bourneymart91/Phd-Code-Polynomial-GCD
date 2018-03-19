function [] = CoefficientPlotting_2Polys
% Plot the coefficients of the polynomials of a set of examples and save
% the figures.


% Create array of example numbers
ex_num_arr = {'18'};





% % Ask user if they wish to save the plots
dialogue_title = 'Save Plots';
opt_1 = 'yes';
opt_2 = 'no';
default_string = 'yes';
bool_savePlots = questdlg('Do you wish to save plots to Figures subfolder?',...
    dialogue_title , opt_1, opt_2, default_string);






% % For each example, plot and save the coefficients of f(x) and g(x)
for i = 1 : 1 : length(ex_num_arr)
    
    ex_num = ex_num_arr{i};
    plotCoefficients(ex_num, bool_savePlots);
    
end



end


function[] = plotCoefficients(ex_num, bool_savePlots)
% Plot coefficients of a given example
%
% % Inputs
%
% ex_num : (String) Example number
%
% bool_savePlots : (String) 'yes' or 'no'


[fx_exact, gx_exact, ~, ~, ~] = Examples_GCD(ex_num);


PlotCoefficients({fx_exact,gx_exact},...
    {'$f(x)$', '$g(x)$'},...
    {'-s', '-o'} ...
    )

switch bool_savePlots
    case 'yes'
        
        file_name = strcat('Example_',ex_num);
        formattype = 'epsc';
        saveas(gcf,file_name,formattype)
        
        formattype = 'png';
        saveas(gcf,file_name,formattype)
        
        formattype = 'jpeg';
        saveas(gcf,file_name,formattype)
        
    case 'no'
        
        
end

end


