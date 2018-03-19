function lambda = GetArithmeticMean(fx, n_k)
% Get arithmetic mean of the entries of the matrix C_{n-k}(f(x))
%
% % Inputs
%
% fx : (Vector) Coefficients of the polynomial f(x)
%
% n_k : (Int) Degree of polynomial v(x)
%
% % Outputs
%
% lambda : (Float) Arithmetic mean of non-zero entries of the partition
% T_{n-k}(f(x)).

global SETTINGS



switch SETTINGS.SYLVESTER_MATRIX_VARIANT 
% Possible Values.
%   * T : 
%   * DT :
%   * DTQ :
%   * TQ : 
%   * DTQ Denominator Removed :
%   * DTQ Rearranged :
    
    
    case 'DTQ'
    
        m = GetDegree(fx);
        
        % Martin's reduced arithmetic mean computation used when computing
        % the arithmetic mean of the entries of a subresultant matrix
        lambda = ((m + n_k + 1) / ((m + 1)^2*(n_k + 1))) * sum(abs(fx));
        
        
    case {'T', 'DT','TQ', 'DTQ Denominator Removed','DTQ Rearranged'}
        
        % Use standard arithmetic mean finding method
        Tf = BuildSubresultant_Partition_2Polys(fx, n_k);
        
        lambda = mean(Tf(Tf~=0));
        
    otherwise 
       
        error('Not a valid Sylvester Matrix build method')
end

end