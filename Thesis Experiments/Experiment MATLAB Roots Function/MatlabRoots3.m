function [] = MatlabRoots3(example_number)
%
% % Inputs
%
% example_number : (Int)

switch example_number
    
    case '1'
        
        exactRoot = 7.123456789;
        max_multiplicity = 40;
        
    case 'Random'
        
        % Random root in interval [a,b]
        a = 0;
        b = 100;
        exactRoot = (b - a).*rand(1,1) + a;
        
        max_multiplicity = 100;
        
    otherwise
end

vAverageDistance = zeros(max_multiplicity,1);

for i = 1 : 1 : max_multiplicity

    multiplicity = i;
    
    % Get vector of exact roots
    vRootsExact = GetRootRepeated(exactRoot, multiplicity);

    % Get computed roots
    calc_roots = roots(poly(vRootsExact));
    
    % Plot computed roots
    figure_name = sprintf('figure_%i',i);
    figure('name',figure_name)
    hold on
    grid on
    scatter(real(calc_roots),imag(calc_roots),'s')
    hold off
    
    
    vDistance = GetDistance(calc_roots, exactRoot);
    vAverageDistance(i) = (mean(vDistance));
    display([calc_roots vDistance])
    
    fx_calc = poly(calc_roots);
    fx_exact = poly(vRootsExact);
    
    vBackwardError(i) = GetError(fx_calc, fx_exact);    
    vForwardError(i) = sum(vDistance);
    
    display(vBackwardError(i));
    display(vForwardError(i));
    
end

sameaxes()

figure()
hold on
plot((vAverageDistance))
hold off

figure()
hold on
plot(vBackwardError,'-s')
hold off

figure()
hold on
plot(vForwardError,'-s')
hold off






end


function vec_root = GetRootRepeated(root, m0)
%
% % Inputs 
%
% a : (Float) Root
%
% m0 : (Int) Multiplicity of root a
%
% % Outputs
%
% vec_a : (Vector) 


% Repeat root m times
vec_root =  root .* ones(1, m0);

end

function [vDistance] = GetDistance(vComputedRoots, exact_root)
% Compute the distance between each computed root and the exact root

% Get number of computed roots
nRoots = length(vComputedRoots);

vDistance = zeros(nRoots,1);

for i = 1:1:nRoots

    % Get calculated root
    calc_root = vComputedRoots(i);
    
    % Get distance between computed root and exact root
    vDistance(i) = norm(exact_root - calc_root);
    
end


end


function [myError] = GetError(fx_comp, fx_exact)

fx_exact_n = fx_exact./norm(fx_exact);

fx_comp_n = fx_comp./norm(fx_comp);

myError = norm(fx_exact_n - fx_comp_n);

end


