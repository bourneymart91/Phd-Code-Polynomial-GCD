function [arr_hxy] = Deconvolve_Bivariate_Batch_Constrained(arr_fxy, vDegt_fxy)
%
%
% % Inputs
%
% arr_fxy : (Array of Matrices) Array of polynomials f_{i}(x,y) in Bernstein form.
%
% vDegt_fxy : (Vector) Vector of total degrees of the polynomials f_{i}(x,y).
%
% % Outputs
%
% arr_hxy : (Array of Matrices) Each cell of the array contains
% coefficients of the polynomial h_{i,j}(x,y)


global SETTINGS

% Get the degree of the polynomials f_{i}(x,y)
vDeg_arr_fxy = vDegt_fxy;

% Get the degree of the polynomials h_{i}(x,y)
vDeg_arr_hxy = diff(vDeg_arr_fxy);

% Get the degree of the polynomials w_{i}(x,y)
vDeg_arr_wxy = diff([vDeg_arr_hxy; 0]);

% Get the multiplicity structure of the polynomial f_{0}(x)
vMult = find(vDeg_arr_wxy~=0);

% Get the number of polynomials in arr_fxy of f_{i}(x,y)
nPolynomials_arr_fxy = size(arr_fxy, 1);

% Get the number of polynomials in the array of h_{i}(x,y)
nPolys_arr_hxy = nPolynomials_arr_fxy - 1;

vDeg_arr_fxy = zeros(nPolynomials_arr_fxy,1);
for i = 1:1:nPolynomials_arr_fxy
    
    vDeg_arr_fxy(i)  = GetDegree_Bivariate(arr_fxy{i});
    
end



% preprocess polynomials f_{i}(x,y)
if SETTINGS.PREPROC_DECONVOLUTIONS
    
    [th1, th2] = GetOptimalTheta(arr_fxy);
    
else
    
    th1 = 1;
    th2 = 1;
    
    
end

% Preprocess polynomials f(x,y) -> f(w_{1}, w_{2})
arr_fww = GetPolynomialArrayWithThetas(arr_fxy, th1, th2);



% %
% %
% Build the LHS Matrix
DT_fwwQ = BuildDTQ_2Polys(arr_fww, vDeg_arr_fxy);

% %
% %
% Build the RHS Vector
rhs_fww = BuildRHS_vec(arr_fww);


x_ls = SolveAx_b(DT_fwwQ, rhs_fww);

nUniquePolys_pxy = unique(vMult);

arr_pww = cell(length(nUniquePolys_pxy), 1);

for i = 1:1:length(nUniquePolys_pxy)
    
    
    mult = nUniquePolys_pxy(i);
    
    % Get degree of p(x,y)
    deg_px = vDegt_fxy(mult) - vDegt_fxy(mult+1);
    
    % Get number of coefficients in p(x,y)
    nCoefficients_px = nchoosek(deg_px + 2, 2);
    nZeros_px = nchoosek(deg_px + 1, 2);
    
    % Get coefficients of p(x,y) from x_ls
    vec_px = x_ls(1 : nCoefficients_px);
    
    % Remove coefficients from solution vector
    x_ls(1 : nCoefficients_px) =[];
    
    
    vec_pxy = ...
        [
        vec_px;
        zeros(nZeros_px, 1)
        ];
    
    arr_pww{i,1} = GetAsMatrix(vec_pxy, deg_px, deg_px);
    
    
    
end


% Get array
arr_hww = Get_hxy(arr_pww, nUniquePolys_pxy);


% Get degree of polynomials h(x,y)
vDeg_arr_hxy = zeros(nPolys_arr_hxy,1);

for i = 1:1:nPolys_arr_hxy
    vDeg_arr_hxy(i) = GetDegree_Bivariate(arr_hww{i});
end



% % Get without thetas
arr_hxy = GetPolynomialArrayWithoutThetas(arr_hww, th1, th2);


end


function LHS_Matrix = BuildDTQ_2Polys(arr_fxy, vDegt_fxy)
%
% % Inputs
%
% arr_fxy : (Array of Matrices)
%
% vDegt_fxy : (Vector) Degree of each polynomial in array
%
% % Outputs
%
% LHS_Matrix : (Matrix)


vDeg_arr_fxy = vDegt_fxy;
vDeg_arr_hxy = diff(vDeg_arr_fxy);
vDeg_wxy = diff([vDeg_arr_hxy; 0]);
vMult = find(vDeg_wxy~=0);


% Get number of distinct polynomials h_{i}(x)
nDistinct_hx = length(vMult);

for i = 1:1:nDistinct_hx
    
    if i>1
        old_mult = vMult(i-1);
    else
        old_mult = 0;
    end
    
    new_mult = vMult(i);
    
    arr_Cf{i} = [];
    
    for j = (old_mult+1+1) : 1 : (new_mult+1)
        
        % Get the degree of previous f(x,y) in the array of polynomials
        deg_fx_prev = vDeg_arr_fxy(j-1);
        
        % Get polynomial f(x,y)
        fxy = arr_fxy{j};
        
        % Get the degree of f(x,y)
        m = vDeg_arr_fxy(j);
        
        % Get the degree of polynomial h_{i}
        deg_hx = deg_fx_prev - m;
        
        % Build the Toeplitz matrix
        arr_Tf{j} = BuildT1(fxy, m, deg_hx);
        
        % Build the diagonal matrix D
        arr_D{j} = BuildD_2Polys(m, deg_hx);
        
        % Stack beneath all other T_f
        arr_Cf{i} = [arr_Cf{i} ; arr_D{j}*arr_Tf{j}];
        
    end
    
    arr_Q{i} = BuildQ1(deg_hx);
    
    arr_DTQ{i} = arr_Cf{i} * arr_Q{i};
end

LHS_Matrix = blkdiag(arr_DTQ{:});

end


function arr_hxy = Get_hxy(arr_pxy, unique_vMult)
%
% % Inputs
%
% arr_pxy : (Array of Matrices)
%
% unique_vMult : (Vector)
%
% % Outputs
%
% arr_hxy :(Array of Matrices)

% Get number of polynomials p_{i}(x,y)
nPolys_arr_pxy = size(arr_pxy,1);

% Initialise a count
count = 1;

for i = 1:1:nPolys_arr_pxy
    
    if i == 1
        nReps = unique_vMult(i);
    else
        nReps = (unique_vMult(i) - unique_vMult(i-1));
    end
    
    for j = 1 : 1 : nReps
        
        arr_hxy{count,1} = arr_pxy{i};
        count = count + 1;
        
    end
    
end
end