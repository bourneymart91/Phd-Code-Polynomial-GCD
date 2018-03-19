function [] = MatlabRoots_StationaryRoots(ex_num)
% Compute the roots of a polynomial, where the multiplicity of one of the
% roots changes.
%
% % Inputs
%
% ex_num : Example Numbers


switch ex_num
    
    case '1'
        vStationary_roots(1) = 1.5;
        vStationary_roots(2) = 4;
        vStationary_roots(3) = 15;
        
        vStationary_root_mult(1) = 7;
        vStationary_root_mult(2) = 10;
        vStationary_root_mult(3) = 10;
        
    case '2'
        vStationary_roots(1) = 1.5;
        
        vStationary_root_mult(1) = 7;
        
    case '3'
        
        vStationary_roots(1) = 2;
        vStationary_roots(2) = 5;
        
        vStationary_root_mult(1) = 7;
        vStationary_root_mult(2) = 2;
end


for i = 0:1:20
    
    myroots = [];
    
    % for each root, repeat it m times
    for j = 1:1:length(vStationary_roots)
        
        % Get root
        my_root = vStationary_roots(j);
        
        % Get root mulitplicity
        my_mult = vStationary_root_mult(j);
        
        % if second root, then vary the multiplicity of the root
        if j == 2
            my_mult = my_mult + i;
        end
        
        % Repeat root m times
        root_vec =  my_root .* ones(1,my_mult);
        
        % Append to vector of roots
        myroots = [myroots root_vec];
    end
    
    
    figure_name = sprintf('Plot-%s',int2str(i));
    
    figure('name',figure_name)
    hold on
    calc_roots = roots(poly(myroots));
    scatter(real(calc_roots),imag(calc_roots),'s')
    
    % Plot exact stationary roots
    for k = 1:1:length(vStationary_roots)
        plot(vStationary_roots(k),0,'*')
    end
    
    hold off
    
    
end


sameaxes()
SaveAllFigures()

end


