function [fx, gx, hx, PolyDx, ux, vx, wx] = ...
    Examples_GCD_FromCoefficients_3Polys(ex_num, ex_num_variant)
%
% % Inputs
%
% ex_num : (String) Example number
%
% ex_num_variant : (String) "a", "b" or "c"
%
% % Outputs
%
% fx : (Vector) Vector of the coefficients of the polynomial f(x)
%
% gx : (Vector) Vector of the coefficients of the polynomial g(x)
%
% hx : (Vector) Vector of the coefficients of the polynomial h(x)
%                   
% dx : (Vector) Vector of the coefficients of the polynomial GCD d(x)
%
% ux : (Vector) Vector of the coefficients of the polynomial u(x)
%
% vx : (Vector) Vector of the coefficients of the polynomial v(x)
%
% wx : (Vector) Vector of the coefficients of the polynomial w(x)

addpath(genpath('../Examples'))


% Get roots and multiplicities from example file
[PolyAx_rm_arr, PolyBx_rm_arr, PolyCx_rm_arr, ...
    PolyDx_rm_arr, PolyPx_rm_arr, PolyQx_rm_arr, Poly_Rx_rm_arr] = ...
    GCD_Examples_Univariate_3Polys(ex_num);



% Get the coefficients of the polynomials f(x), g(x) and h(x) and their GCD
% d(x)
Poly_Ax_sym = GetSymbolicPolyFromSymbolicRoots(PolyAx_rm_arr);
Poly_Bx_sym = GetSymbolicPolyFromSymbolicRoots(PolyBx_rm_arr);
Poly_Cx_sym = GetSymbolicPolyFromSymbolicRoots(PolyCx_rm_arr);
Poly_Dx_sym = GetSymbolicPolyFromSymbolicRoots(PolyDx_rm_arr);


% Get the coefficients of the polynomials f(x), g(x) and d(x) and their GCD
% d(x)
PolyAx = BuildPolyFromRootsSymbolic(PolyAx_rm_arr);
PolyBx = BuildPolyFromRootsSymbolic(PolyBx_rm_arr);
PolyCx = BuildPolyFromRootsSymbolic(PolyCx_rm_arr);
PolyDx = BuildPolyFromRootsSymbolic(PolyDx_rm_arr);


% Get the coefficients of the cofactor polynomials u(x), v(x) and w(x)
PolyPx = BuildPolyFromRootsSymbolic(PolyPx_rm_arr);
PolyQx = BuildPolyFromRootsSymbolic(PolyQx_rm_arr);
PolyRx = BuildPolyFromRootsSymbolic(Poly_Rx_rm_arr);






% % Note - Added extra code to allow poor scaling to be introduced
bad_scaling = false;
if bad_scaling == true
   
    PolyAx = PolyAx * 10^(5);
    PolyBx = PolyBx * 10^(5);
    PolyCx = PolyCx * 10^(-5);
    
end




switch ex_num_variant % Consider the three different polynomial orderings.
    
    case 'a' % Variant a 
        
        fx = PolyAx;
        fx_sym = Poly_Ax_sym;
        
        gx = PolyBx;
        gx_sym = Poly_Bx_sym;
        
        hx = PolyCx;
        hx_sym = Poly_Cx_sym;
        
        ux = PolyPx;
        vx = PolyQx;
        wx = PolyRx;
        
        
    case 'b'
        
        gx = PolyAx;
        gx_sym = Poly_Ax_sym;
        
        fx = PolyBx;
        fx_sym = Poly_Bx_sym;
        hx = PolyCx;
        hx_sym = Poly_Cx_sym;
        vx = PolyPx;
        ux = PolyQx;
        wx = PolyRx;
        
    case 'c'
        
        hx = PolyAx;
        hx_sym = Poly_Ax_sym;
        
        gx = PolyBx;
        gx_sym = Poly_Bx_sym;
        
        fx = PolyCx;
        fx_sym = Poly_Cx_sym;
        
        wx = PolyPx;
        vx = PolyQx;
        ux = PolyRx;
        
    otherwise
        error('Not valid ordering')
end


% Display polynomials f(x) g(x) and h(x)
display(fx_sym)
display(gx_sym)
display(hx_sym)
display(Poly_Dx_sym)


end


