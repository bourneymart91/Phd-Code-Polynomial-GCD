function lambda = GetMean_2Partitions(fx)
% Compute the mean of the non-zero entries of C_{n-k}(f).


global SETTINGS
switch SETTINGS.MEAN_METHOD
    
    case 'Geometric Mean Matlab Method'
        
        % Get geometric mean
        lambda  = geomean(abs(fx(fx~=0)));   
        
    case 'None'
        
        lambda = 1;
        
    otherwise
        
        error('Error SETTINGS.MEAN_METHOD must be valid')
        
end

end