function HYQ = BuildHYQ_SNTLN(dx, m, n, theta)
% Build the matrix HYQ such that H*Y(dx)*Q * [u;v] = [f;g]
%
% % Inputs
%
% dx : (Vector) Coefficients of the polynomial d(x)
%
% m : (Int) Degree of polynomial f(x)
%
% n : (Int) Degree of polynomial g(x)
%
% theta : (Float) Optimal value of theta

HYQ = BuildHYQ1(dx, m, n, theta);

end

function HYQ = BuildHYQ1(dx, m, n, theta)
% Builds the Matrix H.Y.Q where Y_{k} performs a change of variable.
% E(w).d = Y(d).w where d is vector of coefficients of GCD, and w is vector
% of perturbations to f and g
%
% (See Report "APF - General Case - Rearrangement of the matrix vector
% Product")
%
% % Inputs
%
% dx: (Vector) Coefficients of the polynomial d(x)
%
% m : (Int) Degree of polynomial f(x)
%
% n : (Int) Degree of polynomial g(x)
%
% theta : (Float) Optimal value of \theta
%
% % Outputs
%
% HYQ : (Matrix) HYQ


global SETTINGS



if(SETTINGS.BOOL_LOG)
    A = BuildHYPartition_log(dx, m, theta);
    B = BuildHYPartition_log(dx, n, theta);
else
    
    A = BuildHYPartition_nchoosek(dx, m, theta);
    B = BuildHYPartition_nchoosek(dx, n, theta);
    
end

HYQ = blkdiag(A,B);

end

function HY = BuildHYPartition_nchoosek(dx, m, theta)
%
%
% % Inputs.
%
% dx : (Vector) Coefficients of polynomial d(x)
%
% m : (Int) Degree of polynomial f(x)
%
% theta : (Float) Optimal value of \theta
%
% % Outputs
%
% HY : (Matrix)

% Get degree of polynomial d(x)
t = GetDegree(dx);

% Initialise the matrix
HY = zeros(m+1, m-t+1);

% for each column j = 0,\dots,m-t
for j = 0:1:(m-t)
    
    % for each row
    for i = j:1:t+j
        
        
        HY(i+1,j+1) = dx(i-j+1).* theta^(i) ...
            .* nchoosek(i,j) ...]
            .* nchoosek(m-i,t-(i-j)) ;
        
    end
end


% include the denominator
HY = HY./nchoosek(m,t);

end

function A = BuildHYPartition_log(dx, m, theta)
% Build a partition of H*Y
%
% % Inputs
%
%
% dx : (Vector) Coefficients of polynomial d(x)
%
% m : (Int) Degree of polynomial f(x)
%
% theta : (Float) Optimal value of \theta
%
%
% % Outputs
%
% A : (Matrix) Matrix forming a partition of HY
%
%


% Get degree of polynomial d(x)
t = GetDegree(dx);

% Initialise an empty matrix
A = zeros(m+1, m-t+1);

% for each column j of the matrix A
for j = 0:1:(m-t)
    % for each nonzero coefficient.
    for i = j:1:t+j
        
        Num_eval_log = ...
            lnnchoosek(i,j) +...
            lnnchoosek(m-i,t-(i-j));
        
        Num_eval_exp = 10.^Num_eval_log;
        
        A(i+1,j+1) = dx(i-j+1) .* theta^(i) .* Num_eval_exp;
        %A(i+1,j+1) = dx(i-j+1) .* theta^(j) .* Num_eval_exp;
    end
end



Denom_eval_log = lnnchoosek(m,t);
Denom_eval_exp = 10.^Denom_eval_log;
A = A ./ Denom_eval_exp;



end




