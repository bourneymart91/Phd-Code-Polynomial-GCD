function [idx_optColumn] = GetOptimalColumn_Both(S_kk1k2)
% Given two polynomials f(x,y) and g(x,y), and the degree structure of
% their GCD t (t1,t2), pick the optimal column for removal from the 
% Sylvester matrix S_{t_{1},t_{2}}(f,g), where some columns are removed
% from each partition of S_{t1,t2}(f,g) corresponding to the zeros in
% u(x,y) and v(x,y).
%
% Inputs.
%
% Skk1k2 : (Matrix) Sylvester subresultant S_{k,k_{1},k_{2}}(f,g)
%
% Outputs.
%
% idx_optColumn : (Int) Index of optimal column for removal 



% From the given subresultant find the optimal column for removal.
[~,nColumns_S] = size(S_kk1k2);

% Initialise a vector to store residuals by QR decomposition
vResiduals = zeros(nColumns_S,1);

for i = 1:1:nColumns_S
    
    % Get the matrix Ak, which is S_{k} with column c_{k} removed
    Ak = S_kk1k2;
    Ak(:,i) = [];
    
    % Get the column c_{k}
    ck = S_kk1k2(:,i);
    
    % Solve A*x = b
    x_ls = SolveAx_b(Ak,ck);
    
    % Get the residual from this solution
    vResiduals(i) = norm(ck - (Ak*x_ls));
    
end

% Obtain the column for which the residual is minimal.
[~,idx_optColumn] = min(log10(vResiduals));



end