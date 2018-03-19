function vErrors = GetPolynomialArrayErrors(arr_fx_comp, arr_fx_exact)
% Compare each computed h{i} with actual h_{i}
%
% % Inputs
%
% arr_hx_comp : (Array of Vectors) Each vector contains coefficients of the
% polynomial h_{i}(x)
%
% arr_hx_exact : (Array of Vectors) Each vector contains coefficients of
% the polynomial h_{i}(x)

% Get number of polynomials in the array
nPolynomials_arr_hx = size(arr_fx_comp,1);

% Initialise vector to store errors
vErrors = zeros(nPolynomials_arr_hx,1);

%
for i = 1:1:nPolynomials_arr_hx
    
    % Get exact polynomial f_{i}(x)
    fi_exact = arr_fx_exact{i};
    
    % Get computed polynomial f_{i}(x)
    fi_comp = arr_fx_comp{i};
    
    try
    % Get Error
    vErrors(i) = GetPolynomialError(fi_exact, fi_comp);
    
    catch
        vErrors(i) = 1000;
    end
    
    
end


end