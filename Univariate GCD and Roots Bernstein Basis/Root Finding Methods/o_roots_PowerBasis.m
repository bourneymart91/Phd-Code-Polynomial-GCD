function [] = o_roots_Powerbasis(example_number)
% o_roots_powerbasis, takes the polynomial roots, and converts to a
% polynomail in standard power basis, then obtains the roots.

    addpath 'Examples';
    [my_Roots] = Root_Examples(example_number)

    xx = [];
    % For each root
    for i=1:1:length(my_Roots)
       % For each multiplicity
       for j=1:1:my_Roots(i,2)
          xx = [xx my_Roots(i,1)] 
       end
        
    end
    
roots(poly([xx]))



end