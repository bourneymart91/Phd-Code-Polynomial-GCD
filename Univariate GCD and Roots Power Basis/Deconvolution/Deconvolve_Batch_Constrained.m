function [arr_hx] = Deconvolve_Batch_Constrained(arr_fx,vMult)
% Let vMultiplicities be a vector containing the multiplicities of the
% roots of f_{0}(x). vMult = [m_{1}, m_{2} ,..., m_{n}]
% The division f_{0}/f_{1},...,f_{m_{1}-1} / f_{m_{1}} all have the same solution
%
%
% % Inputs.
%
% arr_fx : Array of polynomials f(x)
%
% vMult : Multiplicities of the factors of f_{0}(x) in ascending order.
%
% % Outputs.
%
% arr_hx : Array of polynomials h(x) where h_{i}(x) = f_{i-1}(x) / f_{i}(x)


global SETTINGS

% Get the number of polynomials in the set arr_fx
nPolys_arr_fx = size(arr_fx,1);
nPolys_arr_hx = size(arr_fx,1) - 1;

% Get the degree m_{i} of each of the polynomials f_{i} and store in a
% vector.
vDeg_arr_fx = zeros(nPolys_arr_fx,1);
for i = 1:1:nPolys_arr_fx
    vDeg_arr_fx(i) = GetDegree(arr_fx{i});
end

% Get the degrees n{i} of polynomials h_{i} = f_{i}/f_{i+1}.
vDeg_arr_hx = zeros(nPolys_arr_fx-1,1);
for i = 1:1:nPolys_arr_hx
    vDeg_arr_hx(i) = vDeg_arr_fx(i)-vDeg_arr_fx(i+1);
end

% % Preprocess

if( SETTINGS.PREPROC_DECONVOLUTIONS)    
    theta = GetOptimalTheta(arr_fx,vDeg_arr_fx);
else
    theta = 1;
end


% Get f(\omega) from f(x)
arr_fw = GetPolynomialArrayWithThetas(arr_fx, theta);

% Build LHS Matrix
C_fw = BuildC(arr_fw,vMult);

% % Build the RHS vector
RHS_vec = BuildRHSF(arr_fw);

% % Solve
vec_pw = SolveAx_b(C_fw, RHS_vec);

x_temp = vec_pw;

% Get unique multiplicities of factors of f(x)
unique_vMult = unique(vMult);


arr_pw = cell(length(unique_vMult),1);
for i = 1:1:length(unique_vMult)
    mult = unique_vMult(i);
    deg = GetDegree(arr_fx{mult}) - GetDegree(arr_fx{mult+1});
    
    arr_pw{i} = x_temp(1:deg+1);
    x_temp(1:deg+1) = [];
end


% %
% %
% Get the polynomials p_{i}(x) repeated to give the set of polynomials
% h_{i}(x).

nEntries_arr_px = size(arr_pw, 1);

% Initialise a counter
count = 1;

% Initialise an array to store polynomials h_{i}(\omega)
arr_hw = cell(nPolys_arr_hx,1);

% Get the set of polynomials h_{i}(\omega)
for i = 1:1:nEntries_arr_px
    
    if i == 1
        nReps = unique_vMult(i);
    else
        nReps = (unique_vMult(i) - unique_vMult(i-1));
    end
    
    for j = 1:1:nReps
        arr_hw{count,1} = arr_pw{i,1};
        count = count + 1;
    end
    
end

% Get the set of polynomials h_{i}(x)
arr_hx = GetPolynomialArrayWithoutThetas(arr_hw, theta);



end

function LHS_Matrix = BuildC(arr_fx,vMult)
% %
% %
% Build the LHS Coefficient matrix

% For each distinct hx, build the partition of the Matrix C(f_{m_{i}+1},...,f_{m_{i+1}})
% C(f_{1},...f_{m1}), C(f_{m1+1},...,f_{m2}),...
% C(f_{1},...,f_{m1}) = [T(f1) ; T(f2) ; ... T(fm1)]

% Get number of distinct polynomials in h_{i}
nDistinct_hx = length(vMult);

arr_Cf = cell(1,nDistinct_hx);

for i = 1:1:nDistinct_hx
    
    if i > 1
        old_mult = vMult(i-1);
    else % No previous multiplicity
        old_mult = 0;
    end
    
    new_mult = vMult(i);
    
    arr_Cf{i} = [];
    
    % for each polynomial f_{i} in the interval f_{m_{i-1}+1}...f_{m_{i}}
    
    arr_Tf = cell((new_mult+1) - (old_mult+1+1),1);
    
    for j = (old_mult+1+1) : 1 : (new_mult+1)
        
        % Get the degree of the previous polynomial f_{i-1}(x)
        fx_prev = arr_fx{j-1};
        deg_fx_prev = GetDegree(fx_prev);
        
        % Get the degree of the current polynomial f_{i}(x)
        fx = arr_fx{j};
        deg_fx = GetDegree(fx);
        
        % Get the degree of the polynomial h_{i} = f_{i-1}/f_{i}
        deg_hx = deg_fx_prev - deg_fx;
        
        % Build the Cauchy like matrix T_{m_{i} - m_{i-1}}(f_{i})
        arr_Tf{j} = BuildT1(fx, deg_hx);
        
        % Stack beneath all other T_{f} which are multiplied by [_{i}(x)
        arr_Cf{i} = [arr_Cf{i} ; arr_Tf{j}];
    end
    
    
end

LHS_Matrix = blkdiag(arr_Cf{:});


end