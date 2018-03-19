function [arr_wxy, vDeg_t_wxy] = GetArray_wxy(arr_hxy, vDeg_t_hxy)
%
%
% % Inputs
%
% arr_hxy : Array of polynomials h_{i}(x,y)
%
% vDegt_hxy : Vector of total degrees of h_{i}(x,y)

global SETTINGS

switch SETTINGS.DECONVOLUTION_METHOD_WXY
    
    case 'Separate' % Separate deconvolution
        
        arr_wxy = Deconvolve_Bivariate_Separate(arr_hxy);

        
    case 'Batch'
        
        arr_wxy = Deconvolve_Bivariate_Batch(arr_hxy, vDeg_t_hxy);
        
    case 'Batch With STLN'
        
        arr_wxy = Deconvolve_Bivariate_Batch_With_STLN(arr_hxy, vDeg_t_hxy);
        
    otherwise
        
        display(SETTINGS.DECONVOLUTION_METHOD_WXY)
        
        Options_str = ...
            sprintf([
            '*Separate \n'...
            '*Batch \n'...
            '*Batch With STLN \n'...
            '*Batch Constrained \n'...
            '*Batch Constrained With STLN \n'
            ]);
        error([mfilename ' : ' sprintf('Deconvolution method must be one of the following : \n') Options_str])
        
end

    vDeg_t_wxy = vDeg_t_hxy(1:end-1) - vDeg_t_hxy(2:end);
    
    
    % Set the final w_{x}(i+1) to be equal to h_{x}(i+1)
    arr_wxy{end+1, 1} = arr_hxy{end};
    
    % Get degree
    vDeg_t_wxy(end) = vDeg_t_hxy(end);

for i = 1 : 1 : size(arr_wxy,1)
    fprintf([mfilename ' : ' sprintf('Factors of multiplicity %i :',i) ' \n']);
    factor = arr_wxy{i};
    
    if (length(factor) > 1)
        
        
        fprintf('\n')
        disp(factor)
        %disp(factor./factor(2));
        %disp(factor);
        
    else
        fprintf('\n')
        disp(1)
    end
    
    
end
end