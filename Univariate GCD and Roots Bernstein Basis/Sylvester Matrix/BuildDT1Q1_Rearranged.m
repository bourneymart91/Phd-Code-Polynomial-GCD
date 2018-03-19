function DT1Q1 = BuildDT1Q1_Rearranged(fx,n_k)
% Build the matrix D^{-1}*T_{n-k}(f) * Q_{n-k} in its rearranged format.
%
% Inputs
%
% fx : Coefficients of polynomial f(x,y)
%
% n_k : Degree of polynomial v_{k}(x,y)

global SETTINGS
% If Q is included, use the rearrangment such that each Toeplitz
% matrix has a common divisor in each element.

if(SETTINGS.BOOL_LOG)
    
    % Build Toeplitz Matrix using log version of nchoosek.
    DT1Q1 = BuildDT1Q1_Rearranged_log(fx,n_k);
    
else
    
    % Build Toeplitz Matrix using nchoosek
    DT1Q1 = BuildDT1Q1_Rearranged_nchoosek(fx,n_k);
    
end
end


function DT1Q1 = BuildDT1Q1_Rearranged_log(fx,n_k)
% Build Toeplitz matrix D^{-1}T_{n-k}(f)Q_{n-k}
%
% Inputs.
%
%
% fx :  Coefficients of polynomial f(x) in Bernstein basis.
%
% n_k : Degree of polynomial v(x).

% Get degree of polynomial f(w)
m = GetDegree(fx);

% Build an empty matrix
DT1Q1 = zeros(m+n_k+1,n_k+1);


% For each column of the matrix D*T_{n-k}(f)*Q
for j=0:1:n_k
    
    % for each coefficient in the vector of coefficients of f(x).
    for i = j:1:j+m
        
        % Get the two binomial coefficients in the numerator in terms of
        % logs.
        Numerator_log = lnnchoosek(i,j) + lnnchoosek(m+n_k-i,n_k-j);
        
        % Convert to normal form.
        Numerator_exp = 10.^Numerator_log;
        
        % Multiply binomial evaluation by coefficient f{i-j} from
        % polynomial f.
        DT1Q1(i+1,j+1) =  fx(i-j+1).* Numerator_exp  ;
        
    end
end


% Get log of the binomial coefficient.
Denom_log = lnnchoosek(m+n_k,n_k);

% Get exponent of the log of nchoosek
Denom_exp = 10.^Denom_log;

% Divide each element of T by the common denominator.
DT1Q1 = DT1Q1./Denom_exp;

end

function T = BuildDT1Q1_Rearranged_nchoosek(fx,n_k)
% Build the matrix D^{-1}*T_{n-k}(f)Q_{n-k}
%
% Inputs.
%
% fx :  Coefficients of polynomial f in bernstein basis.
%
% n_k :   Degree of polynomial v(x)


% Get degree of polynomial f(x)
m = GetDegree(fx);

% Build an empty matrix
T = zeros(m+n_k+1,n_k+1);

% for each column j of the Sylvester matrix DTQ
for j=0:1:n_k
    % for each row from i = j,...
    for i = j:1:(j+m)
        T(i+1,j+1) = ...
            fx(i-j+1) ...
            .* nchoosek(i,j) ...
            .* nchoosek(m+n_k-i,n_k-j);
    end
end

scalar = nchoosek(m+n_k,n_k);
T = T./scalar;



end