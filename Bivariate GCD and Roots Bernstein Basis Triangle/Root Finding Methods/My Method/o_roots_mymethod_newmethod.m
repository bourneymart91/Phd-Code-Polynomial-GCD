function [arr_fxy, arr_hxy, arr_wxy, vDeg_t_wxy] = o_roots_mymethod_newmethod(fxy_matrix, M)
% Given a bivariate polynomial compute the roots by differentiating with
% respect to x.
%
% % Inputs
%
%
% fxy_matrix : (Matrix) Coefficients of polynomial f(x,y)
%
% M : Total degree of polynomial f(x,y)
%
% % Outputs
%
% arr_wxy : (Array of Matrices)
%
% vDeg_t_wx : (Vector) Vector containing total degree of each polynomi



% Get the set of polynomials f_{i}(x,y)
[arr_fxy, vDeg_t_arr_fxy] = GetArray_fxy(fxy_matrix, M);

arr_fxy = NormaliseArray(arr_fxy);

% Get array of polynomials h_{i}(x,y)
arr_hxy = GetArray_hxy(arr_fxy, vDeg_t_arr_fxy);

% Get degree of polynomials h(x,y)
vDeg_t_hxy = vDeg_t_arr_fxy(1 : end-1) - vDeg_t_arr_fxy(2 : end);

% Get the set of polynomials w_{i}(x,y)
[arr_wxy, vDeg_t_wxy] = GetArray_wxy(arr_hxy, vDeg_t_hxy);




end


function arr_fxy = NormaliseArray(arr_fxy)

nPolys = length(arr_fxy);


for i = 1 : 1 : nPolys
    
    fxy = arr_fxy{i};
    
    if fxy(1,1) ~= 0
        fxy = fxy ./ fxy(1,1);
    end
    
    arr_fxy{i} = fxy;
    
end

end







