function [ux, vx, wx] = GetQuotients_3Polys(fx, gx, hx, k)
% Given polynomials f(x) and g(x), get the quotient polynomials u(x) and
% v(x) such that f(x)*v(x) = g(x)*u(x).
%
% % Inputs
%
% fx : (Vector) The coefficients of the polynomial f(x)
% 
% gx : (Vector) The coefficients of the polynomial g(x)
%
% hx : (Vector) The coefficients of the polynomial h(x)
% 
% k : (Int) The Degree of common divisor
%
% % Outputs
%
% ux : (Vector) Coefficients of the polynomial u(x)
%
% vx : (Vector) Coefficients of the polynomial v(x)
%
% wx : (Vector) Coefficients of the polynomial w(x)


if (nargin ~= 4)
   erorr('Not enough input arguments') 
end

global SETTINGS

% Get degree of input polynomial g(x) and h(x)
n = GetDegree(gx);
o = GetDegree(hx);

% Build the t^th subresultant
St = BuildSubresultant_3Polys(fx, gx, hx, k);

% Get the optimal column for removal
[~, idx_col] = GetMinDistance(St);

% Remove optimal column
At = St;
At(:, idx_col) = [];

% Get the optimal column c_{t} removed from S_{k}
ct = St(:, idx_col);

% Obtain the solution vector x = [-v;u]
x_ls = SolveAx_b(At, ct);

% Insert a zero into the position corresponding to the index of the optimal
% column so that S(f,g)*vec_x = 0.
vec_x =[
    x_ls(1 : idx_col - 1);
    -1;
    x_ls(idx_col : end);
    ]  ;

% Obtain values for quotient polynomials u and v. still expressed in the
% scaled bernstein basis, including theta.

nCoefficients_vx = n - k + 1;
nCoefficients_wx = o - k + 1;


switch SETTINGS.SYLVESTER_MATRIX_VARIANT
    case {'T'}
        
        vx = vec_x(1 : nCoefficients_vx);
        wx = vec_x(nCoefficients_vx + 1 : nCoefficients_vx + nCoefficients_wx);
        ux = -vec_x(nCoefficients_vx + nCoefficients_wx + 1 : end);
        
    case {'DT'}
        
        
        vx_bi = vec_x(1:nCoefficients_vx);
        wx_bi = vec_x(nCoefficients_vx + 1: nCoefficients_vx + nCoefficients_wx);
        ux_bi = -vec_x(nCoefficients_vx + nCoefficients_wx + 1:end);
        
        ux = GetWithoutBinomials(ux_bi);
        vx = GetWithoutBinomials(vx_bi);
        wx = GetWithoutBinomials(wx_bi);
        
    case {'DTQ'}
        
        vx = vec_x(1 : nCoefficients_vx);
        wx = vec_x(nCoefficients_vx + 1: nCoefficients_vx + nCoefficients_wx);
        ux = -vec_x(nCoefficients_vx + nCoefficients_wx + 1:end);
        
    case {'TQ'}
        
        vx = vec_x(1 : nCoefficients_vx);
        wx = vec_x(nCoefficients_vx + 1 : nCoefficients_vx + nCoefficients_wx);
        ux = -vec_x(nCoefficients_vx + nCoefficients_wx + 1 : end);
        
    case 'DTQ Denominator Removed'
        
        vx = vec_x(1:nCoefficients_vx);
        wx = vec_x(nCoefficients_vx + 1: nCoefficients_vx + nCoefficients_wx);
        ux = -vec_x(nCoefficients_vx + nCoefficients_wx + 1:end);
        
    case 'DTQ Rearranged'
        
        vx = vec_x(1:nCoefficients_vx);
        wx = vec_x(nCoefficients_vx + 1: nCoefficients_vx + nCoefficients_wx);
        ux = -vec_x(nCoefficients_vx + nCoefficients_wx + 1:end);
        
    otherwise 
        error('err')
end



end