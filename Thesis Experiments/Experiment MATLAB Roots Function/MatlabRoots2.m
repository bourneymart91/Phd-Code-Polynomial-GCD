function [] = MatlabRoots2(example_number)
%
% % Inputs
%
% example_number : (Int)


% Root at a with multiplicy m0
% Root at b with multiplicity m1
% Root at c with multiplicity m2

switch example_number
    case 1
        % Location of root 1 and 2
        a = 1;
        b = 2;
        
        % Start position of root 3 (a - s)
        s = 1;
        
        % End position of root 3 (b + t)
        t = 1;

        % Degree of each of the three roots
        m0 = 2;
        m1 = 2;
        m2 = 5;

        % number of positions of root c
        k = 6;
    
    case 2
        
    otherwise
        error('err')
        
end


vec_a = GetRootRepeated(a, m0);
vec_b = GetRootRepeated(b, m1);




for i = 0:1:k
   
    c = (a-s) + ((i * (b + t - a + s)) / k);
    
    vec_c = GetRootRepeated(c, m2);

    vec_allRoots = ...
        [
            vec_a, vec_b, vec_c
        ];
    
    calc_roots = roots(poly(vec_allRoots));
    
    fprintf('Exact Roots \n')
    display(vec_allRoots')
    
    figure_name = sprintf('figure_%i',i);
    figure('name',figure_name)
    hold on
    scatter(real(calc_roots),imag(calc_roots),'s')
    hold off
    
    fprintf('Computed Roots \n')
    display(calc_roots)
    fprintf('\n')
    
end


sameaxes()



    
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