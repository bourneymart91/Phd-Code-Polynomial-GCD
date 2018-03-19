function [idx_optColumn] = GetOptimalColumn_Total(Sk)
% Get index of the optimal column of the Sylvester Matrix S_{k} to be
% removed.
%
% % Inputs.
%
% Sk : (Matrix) Sylvester subresultant matrix S_{k}(f,g)
%
% % Outputs
%
% idx_col : (Int) Index of column to be removed


% From the given subresultant find the optimal column for removal.
[~,nColumns_T1] = size(Sk);

% Take the QR decomposition of the Subresultant
% [Qk,Rk] = qr(Sk);


vResiduals = zeros(nColumns_T1,1);

for i = 1 : 1 : nColumns_T1
    
    %     Sk_temp = Sk;
    %     % Rem
    %     ck = Sk_temp(:,k);
    %     [Q,~] = qrdelete(Qk,Rk,k);
    %     cd = Q'*ck;
    %     d = cd(nCols_T1+1:end,:);
    %     residuals_QR(k) = norm(d);
    
    % Get the matrix Ak, which is S_{k} with column c_{k} removed
    Ak = Sk;
    Ak(:,i) = [];
    
    % Get the column c_{k}
    ck = Sk(:,i);
    
    % Solve A*x = b
    x_ls = SolveAx_b(Ak,ck);
    
    % Get the residual from this solution
    vResiduals(i) = norm(ck - (Ak*x_ls));
    
end

% Obtain the column for which the residual is minimal.
[~,idx_optColumn] = min(log10(vResiduals));


end