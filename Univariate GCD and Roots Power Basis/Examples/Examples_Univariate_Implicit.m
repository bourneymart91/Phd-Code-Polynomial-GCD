function [root_mult_arr_fx] = Examples_Univariate_Implicit(ex_num)
% Given the example number, return the set of roots and corresponding
% multiplicities
%
% Inputs.
%
%
% roots_fx : A matrix where each row consists of [root multiplicity]
% pairs.
%
% Outputs.
%
% root_mult_arr_fx : An array of roots and corresponding multiplicities.

switch ex_num
    case '1'       
        root_mult_arr_fx = ...
            [
            1   1;
            2   2;
            1.7 1;
            ];
    case '2'
        root_mult_arr_fx = ...
            [
            1   1;
            2   2;
            1.5 1;
            ];
        
end

end