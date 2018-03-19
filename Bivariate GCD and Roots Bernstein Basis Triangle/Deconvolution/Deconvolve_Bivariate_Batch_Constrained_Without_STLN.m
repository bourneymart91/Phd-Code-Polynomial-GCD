function [arr_hxy] = Deconvolve_Bivariate_Batch_Constrained_Without_STLN(arr_fxy, vDegt_arr_fxy)
%
%
% % Inputs.
%
% arr_fxy : (Array of Matrices) Array of polynomials f(x,y) in Bernstein form.
%
% vDegt_fxy : (Vector) Vector of degrees of the polynomials f_{i}(x,y).
%
% % Outputs
%
% arr_hxy : (Array of Matrices) Array of polynomials h_{i}(x,y) where each
% cell contains the matrix of coefficients of h_{i}(x,y)

% Get vector of degrees of polynomials f_{i}(x,y)
vDeg_arr_fxy = vDegt_arr_fxy;

% Get vector of degrees of the polynomials h_{i}(x,y)
vDeg_arr_hxy = diff(vDeg_arr_fxy);

% Get vector of degrees of the polynomials w_{i}(x,y)
vDeg_arr_wxy = diff([vDeg_arr_hxy; 0]);

% Get multiplicity structure of the factors of f_{0}(x,y)
vMult = find(vDeg_arr_wxy~=0);


% Get number of polynomials in arr_fxy
nPolys_arr_fxy = size(arr_fxy,1);

vDeg_arr_fxy = zeros(nPolys_arr_fxy,1);
for i = 1:1:nPolys_arr_fxy
    
    vDeg_arr_fxy(i)  = GetDegree_Bivariate(arr_fxy{i});
    
end


% Preprocess
global SETTINGS
if (SETTINGS.PREPROC_DECONVOLUTIONS)
    
    [th1, th2] = GetOptimalTheta(arr_fxy);
    
else
    
    th1 = 1;
    th2 = 1;
    
end

% % Preprocess polynomials f_{i}(x,y)
arr_fww = GetPolynomialArrayWithThetas(arr_fxy,th1,th2);

% % Build the LHS Matrix
DT_fwwQ = BuildDTQ_2Polys(arr_fww, vDeg_arr_fxy);

% %
% %
% Build the RHS Vector
rhs_fww = BuildRHS_vec(arr_fww);


x_ls = SolveAx_b(DT_fwwQ, rhs_fww);

unique_vMult = unique(vMult);

arr_pww = cell(length(unique_vMult),1);

for i = 1:1:length(unique_vMult)
    
    
    mult = unique_vMult(i);
    
    % Get degree of p(x,y)
    deg_px = vDegt_arr_fxy(mult) - vDegt_arr_fxy(mult+1);
    
    % Get number of coefficients in p(x,y)
    nCoefficients_px = nchoosek(deg_px+2, 2);
    
    % Get coefficients of p(x,y) from x_ls
    vec_pww = x_ls(1:nCoefficients_px);
    
    % Remove coefficients
    x_ls(1:nCoefficients_px) =[];
    
    
    nZeros_px = nchoosek(deg_px+1, 2);
    
    vec_pww = ...
        [
        vec_pww;
        zeros(nZeros_px,1)
        ];
    
    arr_pww{i,1} = GetAsMatrix(vec_pww, deg_px, deg_px);
    
    
    
end

% Get number of polynomials 
nPolys_arr_pxy = length(arr_pww);

% Get array of polynomials p_{i}(x,y)
arr_pxy = GetPolynomialArrayWithoutThetas(arr_pww,th1,th2);


count = 1;
for i = 1 : 1 : nPolys_arr_pxy
    
    if i == 1
        nReps = unique_vMult(i);
    else
        nReps = (unique_vMult(i) - unique_vMult(i-1));
    end
    
    for j = 1:1:nReps
        arr_hxy{count,1} = arr_pxy{i};
        count = count + 1;
    end
    
end

end


function LHS_Matrix = BuildDTQ_2Polys(arr_fxy, vDegt_fxy)
%
% % Inputs
%
% arr_fxy : (Array of Matrices)
%
% vDegt_fxy : (Array of Matrices)
%
% % Outputs
%
% LHS_Matrix : (Matrix)

vDeg_f = vDegt_fxy;
vDeg_h = diff(vDeg_f);
vDeg_wxy = diff([vDeg_h; 0]);
vMult = find(vDeg_wxy~=0);
display(vMult);

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
        
        % Get coefficients of previous f(x,y)
        fxy_prev = arr_fxy{j-1};
        
        % Get the degree of previous f(x,y)
        deg_fx_prev = vDeg_f(j-1);
        
        % Get polynomial f(x,y)
        fxy = arr_fxy{j};
        
        % Get the degree of f(x,y)
        deg_fx = vDeg_f(j);
        
        % Get the degree of polynomial h_{i}
        deg_hx = deg_fx_prev - deg_fx;
        
        % Build the cauchy like matrix
        arr_Tf{j} = BuildT1(fxy,deg_fx,deg_hx);
        
        % Build the diagonal matrix D
        arr_D{j} = BuildD(deg_fx,deg_hx);
        
        % Stack beneath all other T_f
        arr_Cf{i} = [arr_Cf{i} ; arr_D{j}*arr_Tf{j}];
        
    end
    
    arr_Q{i} = BuildQ1(deg_hx);
    
    arr_DTQ{i} = arr_Cf{i} * arr_Q{i};
end

LHS_Matrix = blkdiag(arr_DTQ{:});

end
