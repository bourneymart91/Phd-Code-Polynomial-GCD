function [] = CoefficientPlotting_3Polys
% Plot the coefficients of the three polynomials in the three polynomial
% problems found in the example file.




% Initialise an array of examples for which the polynomials coefficients
% are plotted.
ex_num_arr = {'1', '2','3','4','5','6','7','8','9','10','11','12','14','15','16'};





% % Ask user if they wish to save the plots
dialogue_title = 'Save Plots';
opt_1 = 'yes';
opt_2 = 'no';
default_string = 'yes';
bool_savePlots = questdlg('Do you wish to save plots to Figures subfolder?',...
    dialogue_title , opt_1, opt_2, default_string);





% For each example plot the coefficients
for i = 1 : 1 : length(ex_num_arr)

    % Get example number and variant set to 'a'
    ex_num = ex_num_arr{i};
    
    ex_num_variant = 'a';
    
    % Plot coefficients
    plotCoefficients(ex_num, ex_num_variant, bool_savePlots);
    
end



end


function[] = plotCoefficients(ex_num, ex_num_variant, bool_savePlots)
% 
% % Inputs
%
% ex_num : (String) Example number
%
% ex_num_variant : (String) 'a','b' or 'c' Example number variant
% determines the ordering of the polynomials from the example file.
%
% bool_savePlots : (String) 'yes' or 'no'

[fx, gx, hx, ~, ~, ~, ~] = ...
    Examples_GCD_3Polys(ex_num, ex_num_variant);

% Plot coefficients of the polynomials f(x), g(x) and h(x)
PlotCoefficients({fx,gx, hx},...
    {...
    '$\hat{f}(x)$', ...
    '$\hat{g}(x)$', ...
    '$\hat{h}(x)$'...
    },...
    {'-s', '-o','-*'} ...
    )


switch bool_savePlots
    
    case 'yes'
        
        % Set Filename
        file_name = strcat('Example_',ex_num);

        % Save as .eps
        formattype = 'epsc';
        saveas(gcf,file_name,formattype)

        % Save as .fig
        formattype = 'fig';
        saveas(gcf,file_name,formattype)
        
    case 'no'
        
end
end


