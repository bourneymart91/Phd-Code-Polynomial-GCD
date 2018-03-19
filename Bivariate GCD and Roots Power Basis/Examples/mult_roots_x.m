function cellArr = mult_roots_x(root_mult_mat)
% given the root and multiplicity matrix

%                    _______________
% root_mult_mat =   | root  |  mult |
%                   |       |       |
%                   |_______|_______|


count = 1;

cellArr = {};

% for each root in the array
[num_roots,~ ] = size(root_mult_mat);

for i = 1:1:num_roots
    
    % get the multiplicity of the root
    mult = root_mult_mat(i,2);
    
    % get the root
    root = root_mult_mat(i,1);
    
    
    for j = 1:1:mult
        cellArr{count,1} = root_x(root);
        count = count + 1;
    end
    
    
end

end