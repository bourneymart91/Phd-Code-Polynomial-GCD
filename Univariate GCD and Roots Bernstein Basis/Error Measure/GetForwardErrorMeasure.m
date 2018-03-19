function [my_error] = GetForwardErrorMeasure(root_mult_arr_comp, root_mult_arr_exact)
% %
% Used in Root Finding
% Get the error between the set of computed roots and the set of exact
% roots.
%
%
% % Inputs
%
% root_mult_arr_comp : (Matrix)
%
% root_mult_arr_comp : (Matrix)
%
% % Outputs
%
% my_error : (Float)


syms x


nFactors = size(root_mult_arr_exact, 1);

for i = 1 : 1 : nFactors
    
    % Get the factor
    sym_factor = root_mult_arr_exact(i,1);
    
    % Get coefficients in power basis
    try
        pwr_poly = (double((coeffs(sym_factor,x,'All'))))';
    catch
        pwr_poly = (double((coeffs(sym_factor,x))))';
    end
    
    root = - pwr_poly(2) ./ pwr_poly(1);
    
    root_mult_arr_exact(i,1) = root;
    
end


% Sort both matrices based on multiplicity

% Sort the exact matrix
%[values, order] = sort(root_mult_arr_exact(:,2));
%root_mult_arr_exact = root_mult_arr_exact(order,:);
root_mult_arr_exact = sortrows(root_mult_arr_exact,[2,1]);


% Sort the computed matrix
%[values, order] = sort(root_mult_arr_comp(:,2));
%root_mult_arr_comp = root_mult_arr_comp(order,:);
root_mult_arr_comp = sortrows(root_mult_arr_comp,[2,1]);




% Now compare the computed roots
try
    vRoots_exact = double(root_mult_arr_exact(:,1));
    vRoots_comp = double(root_mult_arr_comp(:,1));
    
    vMultiplicity = double(root_mult_arr_exact(:,2));
    
    vDistances = ( vMultiplicity .* abs(vRoots_exact - vRoots_comp));
    
    vDivisors = abs(vRoots_exact);
    
    total_error =  vDistances  ./ vDivisors;
    
    my_error = mean(total_error);
    
catch
    my_error = 1000;
end
end


function [my_error] = GetForwardErrorMeasureArchived(root_mult_arr_comp, root_mult_arr_exact)
% %
% Used in Root Finding
% Get the error between the set of computed roots and the set of exact
% roots.
%
%
% % Inputs
%
% root_mult_arr_comp : (Matrix)
%
% root_mult_arr_comp : (Matrix)
%
% % Outputs
%
% my_error : (Float)


syms x


nFactors = size(root_mult_arr_exact, 1);

for i = 1 : 1 : nFactors
    
    % Get the factor
    sym_factor = root_mult_arr_exact(i,1);
    
    % Get coefficients in power basis
    try
        pwr_poly = double((coeffs(sym_factor,x,'All')))';
    catch
        pwr_poly = double((coeffs(sym_factor,x)))';
    end
    
    root = - pwr_poly(2) ./ pwr_poly(1);
    
    root_mult_arr_exact(i,1) = root;
    
end





% Sort both matrices based on multiplicity

% Sort the exact matrix
%[values, order] = sort(root_mult_arr_exact(:,2));
%root_mult_arr_exact = root_mult_arr_exact(order,:);
root_mult_arr_exact = sortrows(root_mult_arr_exact,[2,1]);


% Sort the computed matrix
%[values, order] = sort(root_mult_arr_comp(:,2));
%root_mult_arr_comp = root_mult_arr_comp(order,:);
root_mult_arr_comp = sortrows(root_mult_arr_comp,[2,1]);


% Check multiplicity structures are equal

vMultExact = root_mult_arr_exact(:,2);
vMultComp = root_mult_arr_comp(:,2);

if sum(abs(vMultExact - vMultComp)) > 0
    fprintf('Multiplicity structures do not match')
    my_error = 1000;
    
else

    % Now compare the computed roots
    try
        vRoots = double(root_mult_arr_exact(:,1));
        vRoots2 = double(root_mult_arr_comp(:,1));
        
        myDistance = abs(vRoots - vRoots2);
        
        total_error = sum( myDistance );
        
        my_error = total_error;
        
    catch
        my_error = 1000;
    end
end
end