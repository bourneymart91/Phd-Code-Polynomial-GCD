function arr_hxy = Deconvolve_Bivariate_Batch(arr_fxy, vDeg_t_fxy, vDeg_x_fxy, vDeg_y_fxy)
% Performs a bivariate deconvolution of polynomials f(x,y) and g(x,y) to
% obtain h(x,y)
%
% % Inputs
%
% arr_fxy : (Array of Matrices) Array of polynomials f_{i}(x,y)
%
% vDeg_t_fxy : (Vector) Total degree of polynomials f_{i}(x,y)
%
% vDeg_x_fxy : (Vector) Degree of polynomials f_{i}(x,y) with respect to x
%
% vDeg_y_fxy : (Vector) Degree of polynomials f_{i}(x,y) with respect to y
%
% % Outputs
%
% arr_hxy : (Array of Matrices) Each array cell contains matrix of 
% coefficients of a polynomial h_{i}(x,y)

global SETTINGS

switch SETTINGS.DEGREE_METHOD
    case 'Total'
        
        arr_hxy = Deconvolve_Bivariate_Batch_Total(arr_fxy, vDeg_t_fxy);
                
    case 'Relative'
        
        arr_hxy = Deconvolve_Bivariate_Batch_Respective(arr_fxy);
        
    case 'Both'
        
        arr_hxy = Deconvolve_Bivariate_Batch_Both(arr_fxy, vDeg_t_fxy, vDeg_x_fxy, vDeg_y_fxy);
        
    otherwise
        
        error('err')
end

end