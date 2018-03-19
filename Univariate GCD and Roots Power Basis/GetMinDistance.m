function [res,index] = GetMinDistance(Sk)
% Given a sylvester matrix get the minimum distance for the removal of a
% column
%
% Inputs.
%
% Sk :

% Get the number of columns in Sk
[~,nColsSk] = size(Sk);

% Initialise a vector to store residual values
residual_vec = zeros(nColsSk,1);

% For each column in S_{k}
for i = 1:1:nColsSk
    
    % Get the matrix Ak, which is S_{k} with column c_{k} removed
    Ak = Sk;
    Ak(:,i) = [];
    
    % Get the column c_{k}
    ck = Sk(:,i);
    
    % Solve A*x = b
    x_ls = SolveAx_b(Ak,ck);
    
    % Get the residual from this solution
    residual_vec(i) = norm(ck - (Ak*x_ls));
    
end

% Get the minimal residual and the index of the corresponding column ck
[res,index] = min(residual_vec);

end