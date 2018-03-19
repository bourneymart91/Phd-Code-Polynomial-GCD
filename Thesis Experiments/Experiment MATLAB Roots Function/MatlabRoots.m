function [] = MatlabRoots(ex_num)
% Compute the roots of the polynomial f(x) = (x-a)^m1(x-b)^m2(x-c)^m3 where
% one of the roots changes, but its multiplicity stays the same.
%
% % Inputs
%
% ex_num : Example Number
%
%
% % Examples
%
% MatlabRoots('1')

switch ex_num
    case '1'
        
        vStationary_roots(1) = 1.5;
        vStationary_roots(2) = 4;
        vStationary_roots(3) = 15;
        
        vStationary_root_mult(1) = 7;
        vStationary_root_mult(2) = 10;
        vStationary_root_mult(3) = 10;

        % Set travelling root start and end point
        travelling_root_start = -1;
        travelling_root_end = 15;
        
        travelling_root_mult = 10;
        
    case '2'
        
        vStationary_roots(1) = 1.5;
       
        vStationary_root_mult(1) = 7;
        
        % Set travelling root start and end point
        travelling_root_start = -1;
        travelling_root_end = 3;
        
        travelling_root_mult = 10;
        
    case '3'
        
        vStationary_roots(1) = 1.5;
        vStationary_roots(2) = 10;
        vStationary_roots(3) = 100;
       
        vStationary_root_mult(1) = 7;
        vStationary_root_mult(2) = 7;
        vStationary_root_mult(3) = 20;
        
        % Set travelling root start and end point
        travelling_root_start = -1;
        travelling_root_end = 110;
        
        travelling_root_mult = 10;
    otherwise 
        error('err')
end


% Get distance travelled by travelling root
travelling_root_distance = travelling_root_end - travelling_root_start;

% Set number of stopping points
nStoppingPoints = 40;

% Get root locations of travelling root
vTravelling_root = zeros(1,nStoppingPoints+1);
for i = 0:1:nStoppingPoints
    vTravelling_root(i+1) = travelling_root_start + ((i/nStoppingPoints) .* travelling_root_distance);
end




% %
% Plot Settings


% % Set limits of y for plotting
% y_lim_lower = -10;
% y_lim_upper = 10;
%
% % Set limits of x for plotting
% x_lim_lower = -10;
% x_lim_upper = 20;



for i = 1:1:length(vTravelling_root)
    
    myroots = [];
    
    % for each root, repeat it m times
    for j = 1:1:length(vStationary_roots)
        
        % Get root
        my_root = vStationary_roots(j);
        
        % Get root mulitplicity
        my_mult = vStationary_root_mult(j);
        
        % Repeat root m times
        root_vec =  my_root .* ones(1,my_mult);
        
        % Append to vector of roots
        myroots = [myroots root_vec];
    end
    
    % Append the moving root
    myroots = [...
        myroots ...
        vTravelling_root(i) .* ones(1,travelling_root_mult)...
        ];
    
    % %
    % Plotting
    
    figure_name = sprintf('Plot-%s',int2str(i));
    
    figure('name',figure_name)
    hold on
    calc_roots = roots(poly(myroots));
    scatter(real(calc_roots),imag(calc_roots),'s')
    
    % Plot exact stationary roots
    for k = 1:1:length(vStationary_roots)
        plot(vStationary_roots(k),0,'*')
    end
    % Plot exact travelling root
    plot(vTravelling_root(i),0,'*');
    hold off
    
    
end

% Ensure all open figures have the same scale
sameaxes()

% Save all figures
savefigs()

end