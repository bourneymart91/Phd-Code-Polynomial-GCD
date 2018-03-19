function H1C1G =  BuildH1C1G(uw,t)
% BuildH1C1G(uw,t)
%
% Build the matrix HCG
%
% Inputs.
%
%
% uw : input polynomial
%
% t : degree of GCD.

global SETTINGS

switch SETTINGS.APF_BUILD_METHOD
    case 'Standard'
        
        m_t = GetDegree(uw);
        m = m_t + t;
        
        % Build the matrix H
        H1 = BuildH1(m);
        
        
        % Build the matrix C1
        C1 = BuildT1(uw,t);
        
        % Build the matrix G
        G  = BuildQ1(t);
        
        % Build the matrix H*C*G
        H1C1G = H1*C1*G;
    case 'Rearranged'
        
        H1C1G = BuildH1C1G_Rearranged(uw,t);
        
    otherwise
        error('SETTINGS.APF_BUILD_METHOD must be either standard or Rearranged')
end


end

function H1C1G = BuildH1C1G_Rearranged(uw,t)
% Build the matrix H1C1G
%
% Inputs.
%
% uw : input polynomial
%
% t : degree of GCD.


% Global Variables
global SETTINGS

if(SETTINGS.BOOL_LOG)
    
    % Use logs
    H1C1G = BuildH1C1G_log(uw,t);
    
elseif (SETTINGS.BOOL_LOG == false)
    
    % Use nchoosek
    H1C1G = BuildH1C1G_nchoosek(uw,t);
    
end

end

function H1C1G = BuildH1C1G_nchoosek(ux, k)
% Build Partition of the HCG matrix using nchoosek
%
%
% Inputs.
%
% ux : (Vector) Coefficients of polynomial u(x)
%
% k : (Int) Degree of GCD
%
%

% Get degree of polynomial u(x)
m_t = GetDegree(ux);

% Get degree of polynomial f(x)
m = m_t + k;

H1C1G = zeros(m+1,k+1);

% for each column 0:1:t
for j = 0:1:k
    %for each row
    for i = j:1:(m_t)+j
        H1C1G(i+1,j+1) = ...
            ux(i-j+1)...
            .* nchoosek(i,j) ...
            .* nchoosek(m-i,k-j);
        
    end
end

% Include Common Denominator in the matrix
H1C1G  = H1C1G ./ nchoosek(m,k);

end

function H1C1G = BuildH1C1G_log(ux, k)
% Build the partition H1C1G where HCG = [H1C1G | H2C2G]
%
% Inputs.
%
%
% uw : (Vector) Coefficients of polynomial u(x)
%
% k : (Int) Index of kth Sylvester subresultant S_{k}(f,g)
%


% Get degree of polynomial uw, deg(u) = m-t.
m_t = GetDegree(ux);

% Get m - the degree of polynomial f.
m = m_t + k;

H1C1G = zeros(m,k);

% for each column 0:1:t
for j = 0:1:k
    %for each row
    for i = j:1:(m_t)+j
        
        Numerator_eval_log = lnnchoosek(i,j) + lnnchoosek(m-i,k-j);
        
        Num_eval_exp = 10.^Numerator_eval_log;
        
        H1C1G(i+1,j+1) = ux(i-j+1) .* Num_eval_exp;
        
        
    end
end


Denom_Eval_log = lnnchoosek(m,k);

Denom_Eval_exp  = 10.^Denom_Eval_log;

H1C1G = H1C1G ./Denom_Eval_exp;


end