function [gm] = GetGeometricMean(fx, n_k)
% This function calculates the geometric mean gm of the terms that contain
% the coefficients of the polynomial c in the kth subresultant matrix
% S(c,d). The integer n is the degree of the polynomial d.
%
% % Inputs
%
% fx : The coefficients of one of the polynomial f(x)
%
% n_k : The degree of the polynomial v(x)
%
%
% % Outputs
%
% gm :  The geometric mean of the terms in the modified Sylvester
%            Sylvester matrix S(f,g)Q.
%



% Dependent on which method is used, 1 - use logs, 0 - use nchoosek
global SETTINGS

switch SETTINGS.BOOL_LOG
    
    case true
        
        gm = GMlog(fx,n_k);
        
    case false
        
        gm = GMnchoosek(fx,n_k);
        
    otherwise
        
        error('BOOL_LOG must be true or false')
end


end


function GM_fx = GMlog(fx, n_k)
% Calculate the Geometric mean of the entries of the Coefficient matrix
% $C_{f}$ which may or may not contain Q. may or may not contain
% denominator.#
%
% % Inputs.
%
% fx : (Vector) Vector containing the coefficients of the polynomial f(x)
%
% n_k : (Int) The degree of polynomial v(x)
%
% % Outputs.
%
% gm : (Float)  Geometric mean of entries of fx in the Syvlester Matrix S_{k}



global SETTINGS


switch SETTINGS.SYLVESTER_MATRIX_VARIANT
    % Values
    %   * T :
    %   * DT :
    %   * DTQ :
    %   * TQ :
    %   * DTQ Denominator Removed :
    %   * DTQ Rearranged :
    
    case 'DTQ'
        % Compute Geometric Mean of D^{-1}T(f,g)Q, which contains three
        % binomial coefficients
        
        % Let.
        %
        % p1 : Product of a_{i}
        %
        % p2 : Product of nchoosek(i,j)
        %
        % p3 : product of nchoosek(m+n-k-(i+j),m-i)
        %
        % p4 : product of nchoosek(m+n-k,m)
        
        % !! Note that p2 = p3
        
        % Calculate the degree of the polynomial.
        m = GetDegree(fx);
        
        % Calculate the absolute value of the coefficients in c.
        fx = abs(fx);
        
        % Calculate part 1 of the geometric mean, the coefficient part.
        % p1_log = (1/(m+1)).*log10(sum(abs(fx)));
        
        p1_log = (1/(m+1)).*sum(log10(abs(fx)));
        
        
        %         % Calculate part 2 of the geometric mean.
        p2_log = 0;
        for j = 0:1:n_k
            for i = 0:1:m
                p2_log = p2_log + lnnchoosek(i+j,j);
            end
        end
        
        p2_log = (1./((n_k+1)*(m+1)))  * (p2_log);
        
        % Calculate part 2 of the geometric mean.
        p3_log = 0;
        for j = 0:1:n_k
            for i = 0:1:m
                p3_log = p3_log + lnnchoosek(m+n_k-(i+j),m-i);
            end
        end
        
        p3_log = (1./((n_k+1)*(m+1)))* (p3_log);
        
        
        
        
        p4_log = lnnchoosek(m+n_k,m);
        
        
        GM_fx = 10.^(p1_log + p3_log + p3_log - p4_log);
        
    case 'DT'
        % Exclude Q from Geometric mean calculations. Split this
        % calculation in to three parts, Numerator_A, the coefficients of
        % fx, Numerator_B : the binomial coefficients corresponding to the
        % a_{i} in the scaled bernstein basis, and Denominator. When
        % Sylvester Matrix is without Q, numerators consist of
        % a_{i}\binom{m}{i}. Denominator changes for each row and is given
        % by \binom{m+n-k}{i+j} where i is the row number and j is the
        % column number, (index from 0).
        
        % Geometric mean is given by $$\prod_{j=0}^{n-k}\prod_{i=0}^{m}$$
        %
        %
        
        % Calculate the degree of the polynomial.
        m = GetDegree(fx);
        
        % Calculate the absolute value of the coefficients in c.
        fx = abs(fx);
        
        % First Dealing with the numerators
        
        Numerator_A = 1;
        Numerator_B = 0;
        
        % For each coefficient a_{i}
        for i = 0:1:m
            Numerator_A = Numerator_A .* fx(i+1);
            Numerator_B = Numerator_B + lnnchoosek(m,i);
        end
        
        GM_Numerator_A = Numerator_A .^(1/(m+1));
        GM_Numerator_B = (1/(m+1)) * Numerator_B;
        GM_Numerator_B = 10.^ GM_Numerator_B;
        GM_Numerator = GM_Numerator_A * GM_Numerator_B;
        
        
        
        Denominator_LOG = 0;
        for i = 0:1:m
            for j = 0:1:n-k
                Denominator_LOG = Denominator_LOG + lnnchoosek(m+n-k,i+j);
            end
        end
        
        
        GM_Denominator_LOG = (1/((m+1)*(n-k+1))) * Denominator_LOG;
        
        GM_Denominator = 10.^GM_Denominator_LOG;
        
        GM_fx = GM_Numerator./GM_Denominator;
        
    case 'T'
        error([mfilename ' : ' 'Code Not Developed']);
        
    case 'TQ'
        error([mfilename ' : ' 'Code Not Developed']);
        
    otherwise
        error('SETTINGS.SYLVESTER_MATRIX_VARIANT not valid');
        
end

end


function GM_fx = GMnchoosek(fx, n_k)
% Get geometric mean of the entries of f(x) in the matrix DTQ using
% nchoosek
%
% Inputs.
%
%
% fx : (Vector) Coefficients of the polynomial f(x)
%
% n_k : (Int) Degree of polynomial v_{k}(x), determines number of columns
% in T_{n-k}(f(x))
%
%
% Outputs.
%
% gm : (Float) Geometric mean of entries of fx in the Syvlester Matrix
% S_{k}
%

% Global Variables.
global SETTINGS


switch SETTINGS.SYLVESTER_MATRIX_VARIANT
    %   * T :
    %   * DT :
    %   * DTQ :
    %   * TQ :
    %   * DTQ Denominator Removed :
    %   * DTQ Rearranged :
    
    
    case 'T'
        
    case 'DT'
        
        % Calculate the degree of the polynomial f.
        m = GetDegree(fx);
        
        fx = abs(fx);
        
        % Get the product of the numerators
        prod_numerator = prod(fx .* GetBinomials(m));
        
        
        % Initialise the product of the denominator at 1.
        prod_denom =1 ;
        
        for i = 0 : 1 : m
            for j = 0 : 1 : n_k
                prod_denom = prod_denom * nchoosek(m+n_k,i+j);
            end
        end
        
        num = prod_numerator .^(1 ./ (m + 1));
        denom = prod_denom .^(1 ./ ((m + 1)*(n_k + 1)));
        
        GM_fx = num/denom;
        
    case 'DTQ'
        
        % Split geometric mean into four parts
        
        % Calculate the degree of the polynomial.
        m = GetDegree(fx);
        p2 = 1;
        
        % Get geometric mean of a_{i} \prod_{i=0}^{m} a_{i} .^(1/m+1) p1 =
        % prod(abs(fx)).^(1./(m+1)) % without using geomean
        
        p1 = geomean(abs(fx));
        
        
        
        % since the product of the binomial coefficient A in the numerator
        % is equal to the product of the binomial coefficient B in the
        % numerator, only calulate this once.
        temp_prod = 1;
        
        for j = 0 : 1 : n_k
            for i = 0 : 1 : m
                temp_prod = temp_prod .* nchoosek(i+j,j);
            end
        end
        
        p2 = temp_prod.^(1./((m+1)*(n_k+1)));
        
        
        % Denominator is included
        p3 = nchoosek(m+n_k,n_k);
        
        
        GM_fx = (p1.*p2.*p2) ./ p3;
        
        
        
    case 'TQ'
        error('err - NOT COMPLETED')
        
    case 'DTQ Denominator Removed'
        error('err - NOT COMPLETED')
        
    case 'DTQ Rearranged'
        error('err - NOT COMPLETED')
        
        
    otherwise
        error('SETTINGS.SYLVESTER_MATRIX_VARIANT is not valid')
end
end

