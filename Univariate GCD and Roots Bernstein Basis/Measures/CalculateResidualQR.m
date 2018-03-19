
function [r] = CalculateResidualQR(ck,Ak)

% This function calculates the residual r of an approximate linear
% algebraic equation whose coefficient matrix is ds and right hand 
% side vector is ck. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Golub & Van Loan Page 263

    [~,n] = size(Ak);
    [Q,~] = qr(Ak);
    cd = Q'*ck;
    d = cd(n+1:end,:);
    r = norm(d);

 
end

