function [fx, gx, hx, dx, ux, vx, wx] = Examples_GCD_FromCoefficients_3Polys(ex_num_var)
%
% % Inputs
%
% ex_num : Example Number
%
% % Outputs
%
% [fx, gx, hx] : Coefficients of the polynomials f(x), g(x) and h(x)
%
% dx : Coefficients of the common divisor d(x)
%
% [ux, vx, wx] : Coefficients of the polynomials u(x), v(x) and w(x)




addpath(genpath('../Examples'))

ex_num = ex_num_var(1 : end - 1);
ex_num_variant = ex_num_var(end);

% Get the rm (Root - Multiplicity) array for each of the polynomials
[PolyA_rm_array, PolyB_rm_array, PolyC_rm_array, ...
    d_rm_array, ...
    Poly1_rm_array, Poly2_rm_array, Poly3_rm_array] = ...
    GCD_Examples_Univariate_3Polys(ex_num);




% Get coefficients of the polynomials f(x), g(x) and h(x)
PolyA = GetCoefficientsFromSymbolicRoots(PolyA_rm_array);
PolyB = GetCoefficientsFromSymbolicRoots(PolyB_rm_array);
PolyC = GetCoefficientsFromSymbolicRoots(PolyC_rm_array);

% Get coefficients of the GCD d(x)
dx = GetCoefficientsFromSymbolicRoots(d_rm_array);

% Get coefficients of u(x), v(x) and w(x)
Poly1 = GetCoefficientsFromSymbolicRoots(Poly1_rm_array);
Poly2 = GetCoefficientsFromSymbolicRoots(Poly2_rm_array);
Poly3 = GetCoefficientsFromSymbolicRoots(Poly3_rm_array);

PolyA_sym = GetSymbolicPolyFromSymbolicRoots(PolyA_rm_array);
PolyB_sym = GetSymbolicPolyFromSymbolicRoots(PolyB_rm_array);
PolyC_sym = GetSymbolicPolyFromSymbolicRoots(PolyC_rm_array);

dx_sym = GetSymbolicPolyFromSymbolicRoots(d_rm_array);




extreme_scaling = false;
if extreme_scaling

    fprintf('NOTE BAD SCALING IS TURNED ON! SEE EXAMPLES_GCD_FROMCOEFFICIENTS_3POLYS.M');
    
    PolyA = PolyA .* 10^(5);
    PolyB = PolyB .* 10^(5);
    PolyC = PolyC .* 10^(-10);

end


% Decide polynomial ordering

switch ex_num_variant
    
    case 'a'
        
        fx = PolyA;
        fx_sym = PolyA_sym;
        gx = PolyB;
        gx_sym = PolyB_sym;
        hx = PolyC;
        hx_sym = PolyC_sym;
        
        ux = Poly1;
        vx = Poly2;
        wx = Poly3;
        
        
    case 'b'
        
        gx = PolyA;
        gx_sym = PolyA_sym;
        fx = PolyB;
        fx_sym = PolyB_sym;
        hx = PolyC;
        hx_sym = PolyC_sym;
        vx = Poly1;
        ux = Poly2;
        wx = Poly3;
        
    case 'c'
        
        hx = PolyA;
        hx_sym = PolyA_sym;
        gx = PolyB;
        gx_sym = PolyB_sym;
        fx = PolyC;
        fx_sym = PolyC_sym;
        wx = Poly1;
        vx = Poly2;
        ux = Poly3;
        
    otherwise
        error('Not valid ordering')
end

% bad_scaling = true;
% 
% if bad_scaling == true
% 
% fx = fx .* 10^(5);
% gx = gx.* 10^(5);
% hx = hx.* 10^(-5);
% 
% end


display(fx_sym)
display(gx_sym)
display(hx_sym)
display(dx_sym)

end

